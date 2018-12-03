
"""
Created on Mon Dec  3 16:26:52 2018

@author: Murat Cem Köse
"""

from setuptools import setup, find_packages

setup(name='SingleR-To-Python',
      version='0.2',
      url='https://github.com/muratcemkose/SingleR-To-Python',
      license='KUL',
      author='Murat Cem Köse',
      author_email='muratckose@gmail.com',
      description='Single cell annotation library',
      packages=find_packages(),
      long_description=open('README.md').read(),
      zip_safe=False)