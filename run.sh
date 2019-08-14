pip install -r requirements.txt

export PYTHONPATH=$PYTHONPATH:.

python autogenes/run.py \
	data/plates \
	2 \
	3 \
	out/ \
	MAON