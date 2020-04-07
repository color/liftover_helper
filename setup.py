import setuptools

setuptools.setup(
    name='LiftoverHelper',
    version='0.1.0',
    author='Anjana Ondov',
    author_email='anju@color.com',
    url='https://github.com/anju24/liftover_helper',
    license='LICENSE.txt',
    description='Pre and post processing of VCF files for liftover.',
    packages=['scripts'],
    install_requires=[
        "pyvcf == 0.6.8",
    ],
    python_requires='>=3.6',
)
