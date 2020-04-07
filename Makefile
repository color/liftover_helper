develop:
	pip install -r requirements.txt

install:
	python setup.py install

clean:
	find tests -name "*.pyc" -exec rm -v {} \;
	find scripts -name "*.pyc" -exec rm -v {} \;
