@rem test.earthmain.gcc.bat: test 32 bit standalone earth.c with main()
@rem
@rem Stephen Milborrow Jan 2008 Durban

@echo test.earthmain.gcc.bat

@mks.cp "C:\bin\R400devdll\i386\R.dll" .
                                @if %errorlevel% neq 0 goto error
@mks.cp "C:\bin\R400devdll\i386\Rblas.dll" .
                                @if %errorlevel% neq 0 goto error
@mks.cp "C:\bin\R400devdll\i386\Riconv.dll" .
                                @if %errorlevel% neq 0 goto error
@mks.cp "C:\bin\R400devdll\i386\Rgraphapp.dll" .
                                @if %errorlevel% neq 0 goto error
@rem you may have to create R.lib and Rblas.lib beforehand
mks.cp "C:\bin\R400devdll\i386\R.lib" .
                                @if %errorlevel% neq 0 goto error
mks.cp "C:\bin\R400devdll\i386\Rblas.lib" .
                                @if %errorlevel% neq 0 goto error

@rem modify the path to include gcc, if needed
@rem only do it if needed
@set | egrep -i "PATH=[^;]*Rtools40.mingw32" >NUL && goto :donesetpath
@echo Modifying path for 32 bit Rtools and R
@set PATH=C:\rtools40\mingw32\bin;^
C:\rtools40\usr\bin;^
C:\Program Files\R\R-4.0.2\bin\i386;^
C:\Program Files\gs\gs9.19\bin;^
%PATH%
:donesetpath

gcc -DSTANDALONE -DMAIN -Wall -pedantic -Wextra -O3 -std=gnu99^
 -I"/a/r/ra/include" -I../../inst/slowtests ../../src/earth.c^
 R.lib Rblas.lib -o earthmain-gcc.exe
                                @if %errorlevel% neq 0 goto error
earthmain-gcc.exe > test.earthmain-gcc.out
                                @rem no errorlevel test, diff will do check for discrepancies
                                @rem @if %errorlevel% neq 0 goto error
@rem we use -w on mks.diff so it treats \r\n the same as \n
mks.diff -w test.earthmain-gcc.out test.earthmain.out.save
                                @if %errorlevel% neq 0 goto error

@rm -f R.dll Rblas.dll Riconv.dll Riconv.dll Rgraphapp.dll R.lib Rblas.lib earthmain-gcc.* test.earthmain-gcc.* *.o
@exit /B 0

:error
@exit /B %errorlevel%
