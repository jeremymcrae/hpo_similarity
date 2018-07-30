"""
Copyright (c) 2015 Wellcome Trust Sanger Institute

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is furnished to do
so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
"""

import io
from setuptools import setup

setup(
    name = "hpo_similarity",
    version = "0.5.1",
    author = "Jeremy McRae",
    author_email = "jeremy.mcrae@sanger.ac.uk",
    description = ("Testing similarity of HPO terms between probands sharing \
        variants in genes."),
    long_description=io.open('README.md', encoding='utf-8').read(),
    long_description_content_type="text/markdown",
    license = "MIT",
    packages=["hpo_similarity"],
    install_requires=['networkx >= 1.8.1',
    ],
    package_data={"hpo_similarity": ['data/hp.obo']},
    classifiers=[
        "Development Status :: 3 - Alpha",
        "License :: OSI Approved :: MIT License",
    ],
    entry_points={'console_scripts': ['hpo_similarity = hpo_similarity.__main__:main']},
    test_suite="tests"
)
