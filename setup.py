from setuptools import setup

__VERSION__ = "0.0.1"


with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="Chemidr",
    version=__VERSION__,
    author="Example Author",
    author_email="author@example.com",
    description="A small example package",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/pypa/sampleproject",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)