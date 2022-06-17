project('ibdools', 'cpp', 
  version : run_command('git', 'describe',  '--tag', check: false).stdout().strip(),
  default_options : ['c_std=c11', 'cpp_std=c++17'])

git_version_h = vcs_tag(command : ['git', 'describe', '--tags'],
            input : 'git_version.txt',
            output :'git_version.h',
            replace_string:'@GIT_VERSION@')

my_inc = include_directories('./include')
my_dep = declare_dependency(include_directories: my_inc)

cxxopts_inc = include_directories('./src/cxxopts/include')
cxxopts_dep = declare_dependency(include_directories: cxxopts_inc)

fmt_dep = dependency('fmt')

htslib_dep = dependency('htslib')

cpps = ['src/common.cpp', 'src/scopedtimer.cpp', 'src/io.cpp', 'src/svsplitter.cpp', 
  'src/tournamenttree.cpp', 'src/ibdfile.cpp', 'src/chromosomes.cpp', 
  'src/genotypes.cpp', 'src/positions.cpp', 'src/metafile.cpp', 'src/ibdsplitters.cpp'
  ]
ibdtools_cpp = ['src/ibdtools.cpp'] + cpps
tests_cpp= ['src/gtest.cpp'] + cpps

# ibdtools
ibdtools_dep = [my_dep, cxxopts_dep, fmt_dep, htslib_dep]
executable('ibdtools', git_version_h, ibdtools_cpp, dependencies: ibdtools_dep)

# tests
gtest_dep = dependency('gtest')
tests_dep = [gtest_dep, my_dep, fmt_dep, htslib_dep]
executable('tests', tests_cpp, dependencies: tests_dep)


