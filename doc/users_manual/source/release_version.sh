#!/bin/sh

old_version="release_version"
release_version="1.1.15"

cat index.rst | sed -e "s/$old_version/$release_version/g" >aaa; mv aaa index.rst 
cat install.rst | sed -e "s/$old_version/$release_version/g" >aaa; mv aaa install.rst 
cat appendix.rst  | sed -e "s/$old_version/$release_version/g" >aaa; mv aaa appendix.rst 
