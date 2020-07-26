from distutils.core import setup
setup(
  name = 'rematbal',
  packages = ['rematbal'],
  version = '0.1',
  license='MIT',
  description = 'Reservoir Engineering Material Balance library',
  author = 'Orkahub Energy',
  author_email = 'orkahub@gmail.com',
  url = 'https://github.com/orkahub/res-eng-mat-bal',
  #download_url = 'https://github.com/user/reponame/archive/v_01.tar.gz',    # I explain this later on
  keywords = ['RESERVOIR', 'ENERGY', 'OIL', 'GAS'],
  install_requires=[
          'pandas',
          'numpy',
          'scipy'
      ],
  classifiers=[
    'Development Status :: 3 - Alpha',
    'Intended Audience :: Developers, Reservoir Engineers',
    'Topic :: Software Development :: Build Tools',
    'License :: OSI Approved :: MIT License',
    'Programming Language :: Python :: 3',
    'Programming Language :: Python :: 3.4',
    'Programming Language :: Python :: 3.5',
    'Programming Language :: Python :: 3.6',
    'Programming Language :: Python :: 3.7',
    'Programming Language :: Python :: 3.8',
  ],
)