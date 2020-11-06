all:
	sudo -H python3 -m pip install -U ortools
	python3 solver.py

install:
	sudo -H python3 -m pip install -U ortools
run:
	python3 solver.py
