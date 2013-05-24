from distutils.core import setup

setup(name="swarm-query",
    version="1.0",
    packages = [''],
    scripts=['swarm_query.py'],
    requires=['pyside', 'numpy', 'matplotlib', 'bsddb3']
    )
