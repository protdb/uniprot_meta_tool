from setuptools import setup, find_packages


def readme():
  with open('README.md', 'r') as f:
    return f.read()


setup(
    name="uniprot_meta_tool",
    description="retrieves and parses metadata from UniProt.org with simple object-based interface",
    long_description=readme(),
    author='Kirill Nikolsky',
    author_email='kirill.s.nikolsky@yandex.ru',
    packages=find_packages(),
    version="0.1.0",
    classifiers=[
        'Programming Language :: Python :: 3.10',
        'Operating System :: OS Independent'
    ],
    install_requires=[
        "requests",
        "jmespath"
    ]
)
