import subprocess
import sys


def executa():
    cmd = ["uv", "run", "python", "cartas_de_controle.py"]
    print("Running:", " ".join(cmd))
    sys.exit(subprocess.run(cmd).returncode)
