import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="hladownload-ramonamerong",
    version="0.0.1",
    author="Ramon van Amerongen",
    author_email="ramonamerong@live.nl",
    description="Python commandline program and module for retrieving HLA alleles, frequencies and eplets.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/pypa/sampleproject",
    project_urls={
        "Bug Tracker": "https://github.com/pypa/sampleproject/issues",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    packages= setuptools.find_packages(include=["hladownload", "hladownload.*"]),
    entry_points = {
        "console_scripts": ['hladownload = hladownload.__main__:main']
    },
    python_requires=">=3.8.5",
    install_requires=[
       'requests>=2.25.1',
       'pandas>=1.2.4',
       'biopython>=1.78',
       'lxml>=4.6.3'
    ]
)