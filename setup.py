from setuptools import setup, find_packages


with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="cenproteo",
    version="0.1.0",
    author="Imiloin, Cannizzaro-reaction, xywawawa",
    description="Some algorithms for finding the essential proteins in a PPI network",
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=find_packages(),
    test_suite="tests",
    python_requires='>=3.8',
    install_requires=[
        "pytest",
    ],
    classifiers=[
        "Development Status :: 3 - Alpha",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
    ],
)