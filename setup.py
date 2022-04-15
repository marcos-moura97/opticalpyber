import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()


setuptools.setup(
    name='opticalpyber',
    version='0.0.1',
    author="Marcos Moura",
    author_email="gabbyru2@gmail.com",
    description="Package with tools to waveguides and optical fibers",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/marcos-moura97/opticalpyber",
    packages=setuptools.find_packages(),
    install_requires=['numpy','scipy'],
    classifiers=[
         "Programming Language :: Python :: 3",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: BSD License",
         "Operating System :: OS Independent",
    ],
 )
