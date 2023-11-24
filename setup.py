from setuptools import setup, find_packages

setup(
    name="Placevent3",
    version="0.1",
    entry_points={
        'console_scripts': [
            'placevent=Placevent3.placevent:main',
        ],
    },
    author="Daniel J. Sindhikara",
    author_email="sindhikara@gmail.com",
    description="Placevent solvent placement software",
    url="https://github.com/mkatouda/Placevent3",
    packages=find_packages(),
    include_package_data=True,
    install_package_data=True,
    install_requires=['numpy'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GPL-3.0 License",
    ],
    python_requires='>=3.8',
)
