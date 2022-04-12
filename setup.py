import setuptools

setuptools.setup(
   name='StochasticSimulations',
   version='0.99',
   description='A collection of tools based on the D-BSSE Stochastic Approaches course',
   author='Nicolò Rossi',
   author_email='olocin.issor@gmail.com',
   install_requires=['wheel', 'numpy', 'scipy'],
   packages=setuptools.find_packages()
)