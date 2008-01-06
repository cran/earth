@rem test.earthmain.gcc.bat: test the standalone earth.c with main()
@rem
@rem Stephen Milborrow Jan 2008 Durban

@echo === test.earthmain.gcc ===
@cp "C:\Program Files\r\R-2.6.1\bin\R.dll" .
@cp "C:\Program Files\r\R-2.6.1\bin\Rblas.dll" .
@cp "C:\Program Files\r\R-2.6.1\bin\iconv.dll" .
@cp "C:\Program Files\r\R-2.6.1\bin\graphapp.dll" .
@cp "C:\Program Files\r\R-2.6.1\bin\Rdll.lib" .
@cp "C:\Program Files\r\R-2.6.1\bin\Rblas.lib" .

@rem modify the path to include gcc, if needed
@set | egrep -i "path=[^;]*Rtools" >NUL && goto :donesetpath
@set path=c:\Rtools\bin;c:\Rtools\MinGW\bin;%PATH%
:donesetpath

@set path=c:\Rtools\bin;c:\Rtools\MinGW\bin;%PATH%
gcc -DSTANDALONE -DMAIN -Wall -pedantic -Wextra -O3 -std=gnu99 -Ic:/a1/r/261/include -I../src/tests ../earth.c Rdll.lib Rblas.lib -o earthmain-gcc.exe
@if %errorlevel% neq 0 goto error:
@earthmain-gcc.exe > test.earthmain-gcc.out
@if %errorlevel% neq 0 goto error:

@echo === Following diff should give no output ===
@rem we use -w on diff so it treats \r\n the same as \n
diff -w test.earthmain.out.save test.earthmain-gcc.out
@if %errorlevel% neq 0 goto error:

@rm -f R.dll Rblas.dll iconv.dll graphapp.dll Rdll.lib Rblas.lib earthmain-gcc.* test.earthmain-gcc.* *.o
@rm -rf Debug
@exit /B 0

:error
@echo error: nmake returned errorlevel %errorlevel%
@exit /B %errorlevel%
