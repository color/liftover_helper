# liftover_helper
Requirements:
Python3

Packages:
pyvcf

# download repo

git clone https://github.com/anju24/liftover_helper.git

cd liftover_helper

python3 -m venv ./venv

# run make
make

# run tests to make sure it works
python -m unittest tests/test_liftover.py
