from setuptools import setup, find_packages

setup(
    name="GeomeTRe",  
    version="0.1.0",         
    description="A package to calculate repeat protein geometrical properties",
    long_description=open("README.md", encoding="utf-8").read(),
    long_description_content_type="text/markdown",
    author="Zarifa Osmanli, Elisa Ferrero",
    author_email="zerifaosmanli@gmail.com",
    url="https://github.com/BioComputingUP/Tandem-repeats-geometry",
    license="GNU General Public License v3.0", 
    packages=find_packages(where="src"),
    package_dir={"": "src"}, 
    include_package_data=True,  
    install_requires=[
        "numpy",
        "pandas",
        "scipy",
        "scikit-learn",
        "biopython",
        "tmtools",
        "scikit-image",
        "requests",
    ],
    entry_points={
        "console_scripts": [
            "geometre=GeomeTRe.main:main",
        ]
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3.0 (GPLv3)",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)


