from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name='ext2TeV',
    version='0.2',
    packages=['ext2TeV', 'ext2TeV.data'],
    package_data={'ext2TeV': ['data/']},
    # packages=find_packages(),
    install_requires=[
        'click',
        'numpy',
        'scipy',
        'matplotlib',
        'astropy',
    ],
    include_package_data=True,
    author="Mireia Nievas",
    author_email="mireia.nievas@iac.es",
    description="A python library to test extrapolation schemes to the VHE regime",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/mireianievas/ext2TeV",
    # packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    # entry_points='''
    # [console_scripts]
    # '''
 )
