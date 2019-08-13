pip install -r requirements.txt

set PYTHONPATH=%PYTHONPATH%;.

python autogenes/run.py data/plates 1 4 out/ MAON

pause