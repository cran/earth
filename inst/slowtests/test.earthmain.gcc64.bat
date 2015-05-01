@rem test.earthmain.gcc.bat: test the standalone earth.c with main()
@rem
@rem TODO I haven't yet been able to get this to work, need 64 bit Rdll.lib and Rblas.lib

@cp "C:/Program Files/R/R-3.2.0/bin/x64/R.dll" .
                                @if %errorlevel% neq 0 goto error
@cp "C:/Program Files/R/R-3.2.0/bin/x64/Rblas.dll" .
                                @if %errorlevel% neq 0 goto error
@cp "C:/Program Files/R/R-3.2.0/bin/x64/Riconv.dll" .
                                @if %errorlevel% neq 0 goto error
@cp "C:/Program Files/R/R-3.2.0/bin/x64/Rgraphapp.dll" .
                                @if %errorlevel% neq 0 goto error
@cp "C:/Program Files/R/R-3.2.0/bin/x64/Rzlib.dll" .
                                @if %errorlevel% neq 0 goto error
@rem you may have to create Rdll.lib and Rblas.lib beforehand
@cp "../../.#/Rdll.lib" .
                                @if %errorlevel% neq 0 goto error
@cp "../../.#/Rblas.lib" .
                                @if %errorlevel% neq 0 goto error
@rem get iconv.dll from /a/r/ra/src/gnuwin32/unicode
@cp "../../.#/Rdll.lib" .
                                @if %errorlevel% neq 0 goto error

@rem modify the path to include gcc, if needed
@set | egrep -i "path=[^;]*Rtools" >NUL && goto donesetpath
@echo Modifying path for Rtools
@set path=D:\Rtools\bin;D:\Rtools\MinGW\bin;%PATH%
:donesetpath

gcc -DSTANDALONE -DMAIN -Wall -pedantic -Wextra -O3 -std=gnu99^
 -m64^
 -I"/a/r/ra/include" -I../../inst/slowtests ../../src/earth.c^
 Rdll.lib Rblas.lib -o earthmain-gcc64.exe
                                @if %errorlevel% neq 0 goto error
@earthmain-gcc.exe > test.earthmain-gcc64.out
                                @if %errorlevel% neq 0 goto error

@rem we use -w on mks.diff so it treats \r\n the same as \n
mks.diff -w test.earthmain-gcc64.out test.earthmain.out64.save
                                @if %errorlevel% neq 0 goto error

@rm -f R.dll Rblas.dll Riconv.dll Riconv.dll Rgraphapp.dll Rzlib.dll Rdll.lib Rblas.lib earthmain-gcc.* test.earthmain-gcc64.* *.o
@exit /B 0

:error
@exit /B %errorlevel%
