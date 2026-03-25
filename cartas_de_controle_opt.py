import os

os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"

import warnings
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import PercentFormatter
from statsmodels.tsa.arima.model import ARIMA
from statsmodels.tsa.arima_process import ArmaProcess

# ============================================================
# Configuração
# ============================================================

SEED = 42

PHI_REAL = 0.5  # valor verdadeiro de Φ1
THETA_REAL = 0.4  # valor verdadeiro de Θ1

tamanhos_amostrais_iniciais = [500]
MC = 100
B = 300

DESVIOS_PHI = [0.0]  # , 0.2]
DESVIOS_THETA = [0.0]  # , 0.2]

SEQ_COLAGENS = np.array([1, 50, 100, 150, 200, 250, 300, 350, 400, 450, 499], dtype=int)

VERBOSE = True

# Ajuste conforme sua máquina.
# Em geral, começar com n_cores - 1 costuma funcionar bem.
N_JOBS = max(1, (os.cpu_count() or 2) - 1)
CHUNKSIZE = 8

MC = 100
B = 300
N_JOBS = 16
CHUNKSIZE = 16


# ============================================================
# Estruturas auxiliares
# ============================================================

@dataclass(slots=True)
class ArmaCoef:
    coef: np.ndarray
    vcov: np.ndarray
    vcov_inv: np.ndarray


# ============================================================
# Funções auxiliares
# ============================================================

# Simula ARMA(1,1)
def simula_arma(n: int, phi: float, theta: float, rng: np.random.Generator, sd: float | None = None) -> np.ndarray:
    """
    Diferença em relação ao R:
    - No R, arima.sim(model = list(ar=phi, ma=theta)) já gera diretamente um ARMA.
    - Em Python/statsmodels, usamos ArmaProcess para simular a série.
    - O statsmodels usa a convenção do polinômio AR com sinal invertido:
      para AR(1), precisamos passar ar = [1, -phi].
    - Para MA(1), usamos ma = [1, theta], compatível com a convenção usual.
    """
    sigma = 1.0 if sd is None else float(sd)

    ar = np.array([1.0, -float(phi)], dtype=float)
    ma = np.array([1.0, float(theta)], dtype=float)

    processo = ArmaProcess(ar, ma)

    # Diferença importante:
    # o gerador interno e detalhes da implementação podem diferir do R,
    # então pequenas diferenças numéricas são esperadas.
    obs = processo.generate_sample(
        nsample=200 + n,
        scale=sigma,
        distrvs=rng.standard_normal,
    )
    return np.asarray(obs[-n:], dtype=float)


# Ajusta ARMA(1,1)
def fit_arma(serie: np.ndarray):
    """
    Diferença em relação ao R:
    - No R, arima(..., order = c(1, 0, 1), include.mean = FALSE, method = "ML")
      devolve um objeto com coef() e vcov().
    - Em Python, usamos statsmodels.tsa.arima.model.ARIMA(...).fit().
    - include.mean = FALSE no R corresponde aqui a trend='n' (sem constante).
    - O método numérico pode não ser idêntico ao R.
    """
    try:
        modelo = ARIMA(serie, order=(1, 0, 1), trend="n")
        return modelo.fit()
    except Exception:
        return None


def extrai_phi_theta(modelo) -> tuple[float | None, float | None]:
    if modelo is None:
        return None, None

    try:
        param_names = list(modelo.param_names)
        params = np.asarray(modelo.params, dtype=float)

        phi_hat = float(params[param_names.index("ar.L1")])
        theta_hat = float(params[param_names.index("ma.L1")])
        return phi_hat, theta_hat
    except Exception:
        return None, None


def arma_coef(modelo) -> ArmaCoef:
    """
    Diferença em relação ao R:
    - coef(modelo) e vcov(modelo) no R viram:
      - modelo.param_names
      - modelo.params
      - modelo.cov_params()
    """
    param_names = list(modelo.param_names)
    params = np.asarray(modelo.params, dtype=float)
    cov = np.asarray(modelo.cov_params(), dtype=float)

    idx_ar1 = param_names.index("ar.L1")
    idx_ma1 = param_names.index("ma.L1")

    coefs = np.array([params[idx_ar1], params[idx_ma1]], dtype=float)
    vcov = cov[np.ix_([idx_ar1, idx_ma1], [idx_ar1, idx_ma1])]
    vcov_inv = np.linalg.inv(vcov)

    return ArmaCoef(coef=coefs, vcov=vcov, vcov_inv=vcov_inv)


def mascara_coef_validos(coefs: np.ndarray) -> np.ndarray:
    phi = coefs[:, 0]
    theta = coefs[:, 1]
    return (
            np.isfinite(phi)
            & np.isfinite(theta)
            & (np.abs(phi) < 0.999)
            & (np.abs(theta) < 0.999)
    )


def amostra_bootstrap_valida(
        mean: np.ndarray,
        cov: np.ndarray,
        B: int,
        rng: np.random.Generator,
        *,
        bloco_inicial: int = 512,
        max_rodadas: int = 20,
) -> np.ndarray:
    """
    Gera B vetores [phi, theta] válidos via normal multivariada, em blocos.

    Diferença importante em relação a uma tradução literal:
    - Em vez de while True por observação, fazemos amostragem vetorial em lotes.
    - Isso é computacionalmente bem mais eficiente.
    """
    aceitos = []
    faltam = B
    bloco = max(bloco_inicial, B)

    for _ in range(max_rodadas):
        amostra = rng.multivariate_normal(mean=mean, cov=cov, size=bloco)
        mask = mascara_coef_validos(amostra)
        validos = amostra[mask]

        if len(validos):
            aceitos.append(validos[:faltam])
            faltam -= min(faltam, len(validos))

        if faltam <= 0:
            break

        bloco = max(bloco, faltam * 3)

    if not aceitos:
        return np.empty((0, 2), dtype=float)

    out = np.vstack(aceitos)
    return out[:B]


def calcula_t2_bootstrap(phi_boot: np.ndarray, theta_boot: np.ndarray):
    boot = np.column_stack([phi_boot, theta_boot])

    medias = boot.mean(axis=0)
    cov = np.cov(boot, rowvar=False, ddof=1)
    cov_inv = np.linalg.inv(cov)

    diffs = boot - medias
    t2 = np.einsum("ij,jk,ik->i", diffs, cov_inv, diffs)

    return medias, cov, cov_inv, t2


# ============================================================
# Worker do bootstrap
# ============================================================

def bootstrap_worker(task: tuple) -> tuple[float, float] | tuple[None, None]:
    """
    Worker isolado para paralelização por processo.

    Entrada:
      (seed, N_INICIAL, sc, serie0, phi_star, theta_star)

    Saída:
      (phi_b, theta_b) ou (None, None)

    Diferença importante:
    - Em paralelo com ProcessPoolExecutor, tudo que vai para o worker
      precisa ser serializável (pickle).
    - Por isso passamos apenas arrays e escalares, não objetos ajustados do statsmodels.
    """
    seed, n_inicial, sc, serie0, phi_star, theta_star = task

    try:
        rng = np.random.default_rng(seed)

        serie_b = simula_arma(n_inicial, float(phi_star), float(theta_star), rng)
        corte = n_inicial - sc
        serie_colada = np.concatenate([serie0[-corte:], serie_b[:sc]])

        fit_b = fit_arma(serie_colada)
        if fit_b is None:
            return None, None

        phi_b, theta_b = extrai_phi_theta(fit_b)
        if phi_b is None or theta_b is None:
            return None, None
        if not np.isfinite(phi_b) or not np.isfinite(theta_b):
            return None, None

        return float(phi_b), float(theta_b)

    except Exception:
        return None, None


# ============================================================
# Execução principal
# ============================================================

def main() -> None:
    rng = np.random.default_rng(SEED)

    if VERBOSE:
        print(f"Iniciando simulações ARMA(1,1) com {MC} iterações Monte Carlo...")
        print(f"Paralelização com {N_JOBS} processos.")

    resultados: list[dict] = []

    # Um pool único ao longo de toda a execução é mais eficiente
    # do que abrir/fechar pool a cada cenário.
    with ProcessPoolExecutor(max_workers=N_JOBS) as executor:
        for mc in range(1, MC + 1):
            if VERBOSE and mc % 10 == 0:
                print(f"Iteração Monte Carlo: {mc}/{MC}")

            for N_INICIAL in tamanhos_amostrais_iniciais:
                # 1) Série base
                serie0 = simula_arma(N_INICIAL, PHI_REAL, THETA_REAL, rng)
                fit0 = fit_arma(serie0)
                if fit0 is None:
                    continue

                try:
                    coef0 = arma_coef(fit0)
                except Exception as e:
                    warnings.warn(f"Falha ao extrair coeficientes da série base: {e}")
                    continue

                seq_colagens_validas = SEQ_COLAGENS[SEQ_COLAGENS < N_INICIAL]

                for desvio_phi in DESVIOS_PHI:
                    for desvio_theta in DESVIOS_THETA:
                        phi_1 = PHI_REAL + desvio_phi
                        theta_1 = THETA_REAL + desvio_theta

                        # Nova série com desvio em φ e θ
                        serie1 = simula_arma(N_INICIAL, phi_1, theta_1, rng)

                        for sc in seq_colagens_validas:
                            if VERBOSE:
                                print(
                                    f"MC={mc}, "
                                    f"N_inicial={N_INICIAL}, "
                                    f"sc={sc}, "
                                    f"desvio_phi={desvio_phi}, "
                                    f"desvio_theta={desvio_theta}"
                                )

                            corte = N_INICIAL - sc

                            # Série colada (controle)
                            serie1_controle = np.concatenate([serie0[-corte:], serie1[:sc]])
                            fit_tmp = fit_arma(serie1_controle)
                            if fit_tmp is None:
                                continue

                            phi_hat, theta_hat = extrai_phi_theta(fit_tmp)
                            if phi_hat is None or theta_hat is None:
                                continue
                            if not np.isfinite(phi_hat) or not np.isfinite(theta_hat):
                                continue

                            # 2) Bootstrap paramétrico: geração dos coeficientes válida no processo principal
                            try:
                                coefs_boot = amostra_bootstrap_valida(
                                    coef0.coef,
                                    coef0.vcov,
                                    B,
                                    rng,
                                    bloco_inicial=max(512, B),
                                    max_rodadas=30,
                                )
                            except Exception as e:
                                warnings.warn(f"Erro na geração de coeficientes bootstrap: {e}")
                                continue

                            if len(coefs_boot) == 0:
                                warnings.warn("Não foi possível gerar coeficientes bootstrap válidos.")
                                continue

                            # Seeds independentes por tarefa
                            task_seeds = rng.integers(
                                low=0,
                                high=np.iinfo(np.uint32).max,
                                size=len(coefs_boot),
                                dtype=np.uint32,
                            )

                            tasks = [
                                (
                                    int(task_seeds[b]),
                                    int(N_INICIAL),
                                    int(sc),
                                    serie0,
                                    float(phi_star),
                                    float(theta_star),
                                )
                                for b, (phi_star, theta_star) in enumerate(coefs_boot)
                            ]

                            # 3) Bootstrap em paralelo
                            results = list(executor.map(bootstrap_worker, tasks, chunksize=CHUNKSIZE))

                            phi_boot_list = []
                            theta_boot_list = []

                            for phi_b, theta_b in results:
                                if phi_b is None or theta_b is None:
                                    continue
                                phi_boot_list.append(phi_b)
                                theta_boot_list.append(theta_b)

                            if len(phi_boot_list) < 5:
                                continue

                            phi_boot = np.asarray(phi_boot_list, dtype=float)
                            theta_boot = np.asarray(theta_boot_list, dtype=float)

                            try:
                                medias_boot, coef_vcov, coef_vcov_inv, t2_boot = calcula_t2_bootstrap(
                                    phi_boot, theta_boot
                                )
                            except Exception as e:
                                warnings.warn(f"Erro ao calcular estatísticas bootstrap: {e}")
                                continue

                            # 4) Limites de controle
                            quant_boot = np.quantile(t2_boot, [0.025, 0.975])

                            # 5) Estatística T² da série de controle
                            diff_controle = np.array([phi_hat, theta_hat], dtype=float) - medias_boot
                            t2_controle = float(diff_controle.T @ coef_vcov_inv @ diff_controle)

                            resultados.append(
                                {
                                    "n0": int(N_INICIAL),
                                    "n1": int(sc),
                                    "desvio_phi": float(desvio_phi),
                                    "desvio_theta": float(desvio_theta),
                                    "fora_de_controle": bool(
                                        (t2_controle < quant_boot[0]) or (t2_controle > quant_boot[1])
                                    ),
                                    "t2.inf": float(quant_boot[0]),
                                    "t2.sup": float(quant_boot[1]),
                                    "t2.controle": t2_controle,
                                }
                            )

    # ============================================================
    # DataFrame final
    # ============================================================

    DADOS_OUT = pd.DataFrame.from_records(
        resultados,
        columns=[
            "n0",
            "n1",
            "desvio_phi",
            "desvio_theta",
            "fora_de_controle",
            "t2.inf",
            "t2.sup",
            "t2.controle",
        ],
    )

    # ============================================================
    # Resumo e gráfico
    # ============================================================

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
        ax.plot(
            df_plot["x_label"].to_numpy(),
            df_plot["proporcao"].to_numpy(),
            marker="o",
            linewidth=1.5,
            label=grupo,
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

    # ============================================================
    # Resumo com IC
    # ============================================================

    DADOS_RESUMO__ = DADOS_RESUMO.copy()
    DADOS_RESUMO__["se"] = np.sqrt(
        DADOS_RESUMO__["proporcao"] * (1.0 - DADOS_RESUMO__["proporcao"]) / MC
    )
    DADOS_RESUMO__["ic_inf"] = np.maximum(
        0.0, DADOS_RESUMO__["proporcao"] - 1.96 * DADOS_RESUMO__["se"]
    )
    DADOS_RESUMO__["ic_sup"] = np.minimum(
        1.0, DADOS_RESUMO__["proporcao"] + 1.96 * DADOS_RESUMO__["se"]
    )

    fig, ax = plt.subplots(figsize=(10, 6))

    for grupo, df_plot in DADOS_RESUMO__.groupby("grupo"):
        df_plot = df_plot.sort_values("x_label")

        x = df_plot["x_label"].to_numpy()
        y = df_plot["proporcao"].to_numpy()
        yerr_inf = (df_plot["proporcao"] - df_plot["ic_inf"]).to_numpy()
        yerr_sup = (df_plot["ic_sup"] - df_plot["proporcao"]).to_numpy()

        ax.plot(x, y, marker="o", linewidth=1.5, label=grupo)
        ax.errorbar(x, y, yerr=[yerr_inf, yerr_sup], fmt="none", capsize=3)

    ax.axhline(0.05, linestyle="dotted")
    ax.set_title("Proporção de Séries Fora de Controle (ARMA(1,1))")
    ax.set_xlabel("Tamanho da Série Colada (n0 + n1)")
    ax.set_ylabel("Proporção Fora de Controle")
    ax.yaxis.set_major_formatter(PercentFormatter(1.0))
    ax.legend(title="Desvios (Φ1; Θ1)", loc="best")
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.show()

    # Export opcional
    # DADOS_OUT.to_csv("py_DADOS_OUT.csv", index=False)
    # DADOS_RESUMO.to_csv("py_DADOS_RESUMO.csv", index=False)


if __name__ == "__main__":
    # No Windows isso é essencial para ProcessPoolExecutor funcionar corretamente.
    main()
