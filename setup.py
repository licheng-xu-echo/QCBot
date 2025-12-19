from setuptools import setup,find_packages
setup(
    name="qcbot",
    version="0.0.4",
    description="Package for a bot for processing quantum calculation data",
    classifiers=["Development Status :: 1 - Pre-Alpha", "Topic :: Scientific/Engineering :: Chemistry"],
    keywords=[],
    url="https://github.com/licheng-xu-echo/QCBot",
    author="Li-Cheng Xu",
    author_email="licheng_xu@zju.edu.cn",
    license="MIT License",
    packages=find_packages(),
    install_package_data=True,
    zip_safe=False,
    install_requires=[],
    package_data={"":["*.csv"]},
)