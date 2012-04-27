@rem test.earthmain.gcc.bat: test the standalone earth.c with main()
@rem
@rem Stephen Milborrow Jan 2008 Durban

@echo === test.earthmain.gcc ===========================================
@cp "C:/Program Files/R/R-2.15.0/bin/i386/R.dll" .
@cp "C:/Program Files/R/R-2.15.0/bin/i386/Rblas.dll" .
@cp "C:/a/r/ra/src/gnuwin32/unicode/iconv.dll" .
@cp "C:/Program Files/R/R-2.15.0/bin/i386/Riconv.dll" .
@cp "C:/Program Files/R/R-2.15.0/bin/i386/Rgraphapp.dll" .
@cp "C:/Program Files/R/R-2.15.0/bin/i386/Rzlib.dll" .
@rem you may have to create Rdll.lib and Rblas.lib beforehand
@cp "../../.#/Rdll.lib" .
@cp "../../.#/Rblas.lib" .

@rem modify the path to include gcc, if needed
@set | egrep -i "path=[^;]*Rtools" >NUL && goto :donesetpath
@set path=c:\Rtools\bin;c:\Rtools\MinGW\bin;%PATH%
:donesetpath

@set path=c:\Rtools\bin;c:\Rtools\MinGW\bin;%PATH%
@gcc -DSTANDALONE -DMAIN -Wall -pedantic -Wextra -O3 -std=gnu99 -I"C:/a/r/ra/include" -I../src/tests ../earth.c Rdll.lib Rblas.lib -o earthmain-gcc.exe
@if %errorlevel% neq 0 goto error
@earthmain-gcc.exe > test.earthmain-gcc.out
@if %errorlevel% neq 0 goto error

@rem we use -w on diff so it treats \r\n the same as \n
diff -w test.earthmain-gcc.out test.earthmain.out.save
@if %errorlevel% neq 0 goto error

@rm -f R.dll Rblas.dll Riconv.dll Riconv.dll Rgraphapp.dll Rzlib.dll Rdll.lib Rblas.lib earthmain-gcc.* test.earthmain-gcc.* *.o
@exit /B 0

:error
@echo error: errorlevel %errorlevel%
@exit /B %errorlevel%
