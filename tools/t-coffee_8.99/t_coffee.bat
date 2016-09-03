@echo off

set ROOTDIR="%~d0%~p0"
set DISK=%~d0

%DISK%
cd %ROOTDIR%
bin\sh bin/t_coffee %*
