"""Run the analysis and figures using independently obtained authorized inputs."""
from pathlib import Path
import os,subprocess,sys

def main():
    if not os.environ.get('SEER_EXPORT_DIR'):
        raise SystemExit('Set SEER_EXPORT_DIR to the authorized local export directory.')
    root=Path(__file__).resolve().parent
    for script in ['analyze.py','make_figures.py']:
        subprocess.run([sys.executable,str(root/script)],check=True)

if __name__=='__main__':
    main()
