set -e

SCRIPT_DIR="$(dirname "$0")"
if [ ! -d ~/.venv ]; then
	python3.12 -m venv ~/.venv
fi
source ~/.venv/bin/activate
pip3 install --no-cache-dir -r "$(dirname "$SCRIPT_DIR")"/requirements.txt
deactivate
