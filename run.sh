pip install -r requirements.txt

export PYTHONPATH=$PYTHONPATH:.

python autogenes/run.py \
	data/plates \
	1 \
	4 \
	out/ \
	MAON