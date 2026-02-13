import subprocess
from pathlib import Path

examples_dir = Path("examples")

for script in examples_dir.glob("*.py"):
    print(f"Running {script.name}")

    subprocess.run(
        ["python", str(script)],
        check=True
    )