from setuptools import setup, find_packages

setup(
    name="geometre",
    version="0.1.0",
    description="Calculate repeat protein geometrical properties",
    long_description=open("README.md", encoding="utf-8").read(),
    long_description_content_type="text/markdown",
    author="Zarifa Osmanli, Elisa Ferrero, Damiano Piovesan",
    author_email="zerifaosmanli@gmail.com, damiano.piovesan@unipd.it",
    url="https://github.com/BioComputingUP/GeomeTRe",
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
    extras_require = {
        'draw':  ["pymol"]
    },
    entry_points={
        "console_scripts": [
            "geometre=geometre.main:main",
        ]
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3.0 (GPLv3)",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)


