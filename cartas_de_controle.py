import warnings

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import PercentFormatter
from scipy.stats import multivariate_normal
from statsmodels.tsa.arima.model import ARIMA, ARIMAResults

np.random.seed(42)

# ------------------------------
# Parâmetros gerais
# ------------------------------
PHI_REAL = 0.5  # valor verdadeiro de Φ1
THETA_REAL = 0.4  # valor verdadeiro de Θ1

tamanhos_amostrais_iniciais = [500]
MC = 100
B = 300

DESVIOS_PHI = [0.0, 0.2]  # , 0.6)
DESVIOS_THETA = [0.0, 0.2]  # , 0.6)

# passo = 20
# bc = [1, 5]
# SEQ_COLAGENS = [y for x in range(0, 101, passo) for y in (bc[0] + x, bc[1] + x)]
# SEQ_COLAGENS = SEQ_COLAGENS + [135, 150, 165, 175, 185, 199]
SEQ_COLAGENS = [1, 50, 100, 150, 200, 250, 300, 350, 400, 450, 499]

DADOS_OUT = pd.DataFrame(
    columns=[
        "n0",
        "n1",
        "desvio_phi",
        "desvio_theta",
        "fora_de_controle",
        "t2.inf",
        "t2.sup",
        "t2.controle",
    ]
)


# ------------------------------
# Funções auxiliares
# ------------------------------

# Simula ARMA(1,1)
def simula_arma(n: int, phi: float, theta: float, sd: float | None = None) -> np.ndarray:
    """
    Diferença em relação ao R:
    - No R, arima.sim(model = list(ar=phi, ma=theta)) já gera diretamente um ARMA.
    - Em Python/statsmodels, usamos ArmaProcess para simular a série.
    - O statsmodels usa a convenção do polinômio AR com sinal invertido:
      para AR(1), precisamos passar ar = [1, -phi].
    - Para MA(1), usamos ma = [1, theta], compatível com a convenção usual.
    """
    from statsmodels.tsa.arima_process import ArmaProcess

    sigma = 1.0 if sd is None else float(sd)

    ar = np.array([1.0, -float(phi)])
    ma = np.array([1.0, float(theta)])

    processo = ArmaProcess(ar, ma)
    obs = processo.generate_sample(nsample=200 + n, scale=sigma)
    return obs[-n:]


# Ajusta ARMA(1,1)
def fit_arma(serie: np.ndarray) -> ARIMAResults | None:
    """
    Diferença em relação ao R:
    - No R, arima(..., order = c(1, 0, 1), include.mean = FALSE, method = "ML")
      devolve um objeto com coef() e vcov().
    - Em Python, usamos statsmodels.tsa.arima.model.ARIMA(...).fit().
    - include.mean = FALSE no R corresponde aqui a trend='n' (sem constante).
    - O método exato de estimação interna pode não ser idêntico ao R, então
      pequenas diferenças numéricas são esperadas mesmo com a mesma semente.
    """
    try:
        modelo = ARIMA(serie, order=(1, 0, 1), trend="n")
        resultado = modelo.fit()
        return resultado
    except Exception as e:
        warnings.warn(f"Erro ao ajustar ARMA(1,1): {e}")
        return None


# Diferença importante em relação à checagem por raízes:
# para ARMA(1,1), sob a parametrização usada aqui, as condições
# práticas são |phi| < 1 (estacionariedade) e |theta| < 1 (invertibilidade).
# Embora em teoria isso possa ser expresso via raízes dos polinômios,
# usar diretamente essas desigualdades é mais simples e evita erro de convenção
# entre a forma do polinômio em R/statsmodels e a forma implementada manualmente.

# Testa estacionaridade (para φ)
def ar_valido(phi: float) -> bool:
    # AR(1) estacionário <=> |phi| < 1
    return np.isfinite(phi) and abs(phi) < 1.0


# Testa invertibilidade (para θ)
def ma_valido(theta: float) -> bool:
    # MA(1) invertível <=> |theta| < 1
    return np.isfinite(theta) and abs(theta) < 1.0


# Extrai coeficientes e matriz de covariância
def arma_coef(modelo: ARIMAResults) -> dict[str, np.ndarray]:
    """
    Diferença em relação ao R:
    - coef(modelo) e vcov(modelo) no R viram, em statsmodels:
      - modelo.param_names
      - modelo.params
      - modelo.cov_params()
    - Como o vetor de parâmetros pode conter nomes diferentes dependendo da versão,
      fazemos o mapeamento por nome.
    """
    param_names = list(modelo.param_names)
    params = np.asarray(modelo.params, dtype=float)
    cov = np.asarray(modelo.cov_params(), dtype=float)

    idx_ar1 = param_names.index("ar.L1")
    idx_ma1 = param_names.index("ma.L1")

    coefs = np.array([params[idx_ar1], params[idx_ma1]], dtype=float)
    vcov = cov[np.ix_([idx_ar1, idx_ma1], [idx_ar1, idx_ma1])]

    return {
        "coef": coefs,
        "vcov": vcov,
        "vcov_inv": np.linalg.inv(vcov),
    }


def extrai_phi_theta(modelo: ARIMAResults | None) -> tuple[float | None, float | None]:
    """
    Extrai os coeficientes AR(1) e MA(1) do resultado ajustado.
    """
    if modelo is None:
        return None, None

    try:
        param_names = list(modelo.param_names)
        params = np.asarray(modelo.params, dtype=float)

        phi_hat_ = float(params[param_names.index("ar.L1")])
        theta_hat_ = float(params[param_names.index("ma.L1")])
        return phi_hat_, theta_hat_
    except Exception as e:
        warnings.warn(f"Erro ao extrair coeficientes ARMA: {e}")
        return None, None


# ------------------------------
# Loop Monte Carlo
# ------------------------------

print(f"Iniciando simulações ARMA(1,1) com {MC} iterações Monte Carlo...")

for mc in range(1, MC + 1):
    if mc % 10 == 0:
        print(f"Iteração Monte Carlo: {mc}/{MC}")

    for N_INICIAL in tamanhos_amostrais_iniciais:
        # 1) Série base
        serie0 = simula_arma(N_INICIAL, PHI_REAL, THETA_REAL)
        fit0 = fit_arma(serie0)
        if fit0 is None:
            continue

        coef0 = arma_coef(fit0)

        for desvio_phi in DESVIOS_PHI:
            for desvio_theta in DESVIOS_THETA:
                # Nova série com desvio em φ e θ
                serie1 = simula_arma(
                    N_INICIAL,
                    PHI_REAL + desvio_phi,
                    THETA_REAL + desvio_theta,
                )

                for sc in SEQ_COLAGENS:
                    if sc >= N_INICIAL:
                        break

                    print(
                        f"MC={mc}, "
                        f"N_inicial={N_INICIAL}, "
                        f"sc={sc}, "
                        f"desvio_phi={desvio_phi}, "
                        f"desvio_theta={desvio_theta}"
                    )

                    # Série colada (controle)
                    serie1_controle = np.concatenate(
                        [serie0[-(N_INICIAL - sc):], serie1[:sc]]
                    )
                    assert len(serie1_controle) == N_INICIAL

                    fit_tmp = fit_arma(serie1_controle)
                    if fit_tmp is None:
                        continue

                    phi_hat, theta_hat = extrai_phi_theta(fit_tmp)
                    if phi_hat is None or theta_hat is None:
                        continue
                    if np.isnan(phi_hat) or np.isnan(theta_hat):
                        continue

                    # 2) Bootstrap paramétrico
                    phi_boot = np.full(B, np.nan, dtype=float)
                    theta_boot = np.full(B, np.nan, dtype=float)

                    for b in range(B):
                        try:
                            # 'Sigma' is not positive definite
                            while True:
                                coef_star = multivariate_normal.rvs(
                                    mean=coef0["coef"],
                                    cov=coef0["vcov"],
                                )
                                coef_star = np.asarray(coef_star, dtype=float)

                                phi_star = float(coef_star[0])
                                theta_star = float(coef_star[1])

                                if ar_valido(phi_star) and ma_valido(theta_star):
                                    break
                        except Exception as e:
                            warnings.warn(f"Erro na geração de coeficientes bootstrap: {e}")
                            continue

                        serie_b = simula_arma(N_INICIAL, phi_star, theta_star)
                        serie_colada = np.concatenate(
                            [serie0[-(N_INICIAL - sc):], serie_b[:sc]]
                        )
                        assert len(serie_colada) == N_INICIAL

                        fit_b = fit_arma(serie_colada)
                        if fit_b is None:
                            phi_boot[b] = np.nan
                            theta_boot[b] = np.nan
                            continue

                        phi_b, theta_b = extrai_phi_theta(fit_b)

                        if phi_b is None or theta_b is None:
                            print("Ajuste ARMA falhou no bootstrap, pulando iteração.")
                            phi_boot[b] = np.nan
                            theta_boot[b] = np.nan
                            continue

                        phi_boot[b] = phi_b
                        theta_boot[b] = theta_b

                    phi_boot = phi_boot[~np.isnan(phi_boot)]
                    theta_boot = theta_boot[~np.isnan(theta_boot)]

                    if len(phi_boot) < 5 or len(theta_boot) < 5:
                        continue

                    coef_vcov = np.array(
                        [
                            [np.var(phi_boot, ddof=1), np.cov(phi_boot, theta_boot, ddof=1)[0, 1]],
                            [np.cov(theta_boot, phi_boot, ddof=1)[0, 1], np.var(theta_boot, ddof=1)],
                        ],
                        dtype=float,
                    )

                    coef_vcov_inv = np.linalg.inv(coef_vcov)
                    phi_m = float(np.mean(phi_boot))
                    theta_m = float(np.mean(theta_boot))

                    t2 = np.empty(len(phi_boot), dtype=float)
                    for i in range(len(phi_boot)):
                        diff_b = np.array([phi_boot[i], theta_boot[i]], dtype=float) - np.array(
                            [phi_m, theta_m], dtype=float
                        )
                        t2[i] = float(diff_b.T @ coef_vcov_inv @ diff_b)

                    # 3) Limites de controle
                    quant_boot = np.quantile(t2, [0.025, 0.975])

                    # 4) Estatística T²
                    diff_controle = np.array([phi_hat, theta_hat], dtype=float) - np.array(
                        [phi_m, theta_m], dtype=float
                    )
                    t2_controle = float(diff_controle.T @ coef_vcov_inv @ diff_controle)

                    DADOS_OUT.loc[len(DADOS_OUT)] = {
                        "n0": N_INICIAL,
                        "n1": sc,
                        "desvio_phi": desvio_phi,
                        "desvio_theta": desvio_theta,
                        "fora_de_controle": (
                                t2_controle < quant_boot[0] or t2_controle > quant_boot[1]
                        ),
                        "t2.inf": float(quant_boot[0]),
                        "t2.sup": float(quant_boot[1]),
                        "t2.controle": t2_controle,
                    }

# ------------------------------
# Resumo e gráfico
# ------------------------------

DADOS_RESUMO = (
    DADOS_OUT.groupby(["n0", "n1", "desvio_phi", "desvio_theta"], as_index=False)
    .agg(proporcao=("fora_de_controle", "mean"))
)
DADOS_RESUMO["x_label"] = DADOS_RESUMO["n1"] + DADOS_RESUMO["n0"]
DADOS_RESUMO["grupo"] = (
        DADOS_RESUMO["desvio_phi"].astype(str) + ";" + DADOS_RESUMO["desvio_theta"].astype(str)
)

fig, ax = plt.subplots(figsize=(10, 6))

for grupo, df_plot in DADOS_RESUMO.groupby("grupo"):
    df_plot = df_plot.sort_values("x_label")
    ax.plot(df_plot["x_label"], df_plot["proporcao"], marker="o", linewidth=1.5, label=grupo)

ax.axhline(0.05, linestyle="dotted")
ax.set_title("Proporção de Séries Fora de Controle (ARMA(1,1))")
ax.set_xlabel("Tamanho da Série Colada (n0 + n1)")
ax.set_ylabel("Proporção Fora de Controle")
ax.yaxis.set_major_formatter(PercentFormatter(1.0))
ax.legend(title="Desvios (Φ1; Θ1)", loc="best")
plt.xticks(rotation=45)
plt.tight_layout()
plt.show()

DADOS_RESUMO__ = (
    DADOS_OUT.groupby(["n0", "n1", "desvio_phi", "desvio_theta"], as_index=False)
    .agg(proporcao=("fora_de_controle", "mean"))
)
DADOS_RESUMO__["x_label"] = DADOS_RESUMO__["n0"] + DADOS_RESUMO__["n1"]
DADOS_RESUMO__["se"] = np.sqrt(
    DADOS_RESUMO__["proporcao"] * (1.0 - DADOS_RESUMO__["proporcao"]) / MC
)
DADOS_RESUMO__["ic_inf"] = np.maximum(0.0, DADOS_RESUMO__["proporcao"] - 1.96 * DADOS_RESUMO__["se"])
DADOS_RESUMO__["ic_sup"] = np.minimum(1.0, DADOS_RESUMO__["proporcao"] + 1.96 * DADOS_RESUMO__["se"])
DADOS_RESUMO__["grupo"] = (
        DADOS_RESUMO__["desvio_phi"].astype(str) + ";" + DADOS_RESUMO__["desvio_theta"].astype(str)
)

fig, ax = plt.subplots(figsize=(10, 6))

for grupo, df_plot in DADOS_RESUMO__.groupby("grupo"):
    df_plot = df_plot.sort_values("x_label")
    ax.plot(df_plot["x_label"], df_plot["proporcao"], marker="o", linewidth=1.5, label=grupo)
    ax.errorbar(
        df_plot["x_label"],
        df_plot["proporcao"],
        yerr=[
            df_plot["proporcao"] - df_plot["ic_inf"],
            df_plot["ic_sup"] - df_plot["proporcao"],
        ],
        fmt="none",
        capsize=3,
    )

ax.axhline(0.05, linestyle="dotted")
ax.set_title("Proporção de Séries Fora de Controle (ARMA(1,1))")
ax.set_xlabel("Tamanho da Série Colada (n0 + n1)")
ax.set_ylabel("Proporção Fora de Controle")
ax.yaxis.set_major_formatter(PercentFormatter(1.0))
ax.legend(title="Desvios (Φ1; Θ1)", loc="best")
plt.xticks(rotation=45)
plt.tight_layout()
plt.show()

if __name__ == "__main__":
    print("Script finalizado com sucesso!")
