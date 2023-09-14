from setuptools import setup, find_namespace_packages


setup(
    name='core-typer',
    version='0.1.0-alpha',
    packages=find_namespace_packages(),
    entry_points={
        "console_scripts": [
            "core-typer = core_typer.__main__:main",
        ]
    },
    scripts=[],
    package_data={
    },
    install_requires=[
    ],
    description='A cgMLST Typing Tool',
    url='https://github.com/dfornika/core-typer',
    author='Dan Fornika',
    author_email='dfornika@gmail.com',
    include_package_data=True,
    keywords=[],
    zip_safe=False
)
