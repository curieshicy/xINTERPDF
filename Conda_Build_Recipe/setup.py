from setuptools import setup

requirements = [
    # package requirements go here
]

setup(
    name='xinterpdf',
    version='0.1.0',
    description="a GUI for analyzing X-ray PDF data for organics",
    author="Chenyang Shi",
    author_email='chenyang.shi@abbvie.com',
    url='https://github.com/curieshicy/xINTERPDF',
    packages=['xinterpdf'],
    include_package_data=True,
    package_data={'xinterpdf': ['*.gif']},
    entry_points={
        'console_scripts': [
            'xinterpdf=xinterpdf.cli:main'
        ]
    },
    install_requires=requirements,
    keywords='xinterpdf',
    classifiers=[
        'Programming Language :: Python :: 2.7'    ]
)
