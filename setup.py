from setuptools import setup

setup(name='clumo',
      version='0.0',
      description='Mock galaxy cluster generation tool',
      url='https://github.com/lavanyanemani1996/clumo.git',
      install_requires=['numpy', 'math', 'astropy', 'hmf', 'scipy','random',
			'time', 'multiprocessing', 'itertools','xspec'],
      packages=['myclumo'],
      zip_safe=False)
