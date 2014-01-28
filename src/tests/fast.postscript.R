# fast.postscript.R:
#
# Function to open a postscript file that does not use the sRGB handling
# introduced in R version 2.13-0.  The resulting postscript file will display
# several times faster than a file that uses the sRGB handling.

fast.postscript <- function(file= ifelse(onefile, "Rplots.ps", "Rplot%03d.ps"),
                            onefile=TRUE, ...)
{
    prolog <- c(
    "/gs  { gsave } def",
    "/gr  { grestore } def",
    "/ep  { showpage gr gr } def",
    "/m   { moveto } def",
    "/l  { rlineto } def",
    "/np  { newpath } def",
    "/cp  { closepath } def",
    "/f   { fill } def",
    "/o   { stroke } def",
    "/c   { newpath 0 360 arc } def",
    "/r   { 4 2 roll moveto 1 copy 3 -1 roll exch 0 exch rlineto 0 rlineto -1 mul 0 exch rlineto closepath } def",
    "/p1  { stroke } def",
    "/p2  { gsave bg fill grestore newpath } def",
    "/p3  { gsave bg fill grestore stroke } def",
    "/p6  { gsave bg eofill grestore newpath } def",
    "/p7  { gsave bg eofill grestore stroke } def",
    "/t   { 5 -2 roll moveto gsave rotate",
    "       1 index stringwidth pop",
    "       mul neg 0 rmoveto show grestore } def",
    "/ta  { 4 -2 roll moveto gsave rotate show } def",
    "/tb  { 2 -1 roll 0 rmoveto show } def",
    "/cl  { grestore gsave newpath 3 index 3 index moveto 1 index",
    "       4 -1 roll lineto  exch 1 index lineto lineto",
    "       closepath clip newpath } def",

    ### original R 2.13-0 code
    # "/srgb { [ /CIEBasedABC",
    # "          << /DecodeLMN",
    # "               [ { dup 0.03928 le",
    # "                        {12.92321 div}",
    # "                        {0.055 add 1.055 div 2.4 exp }",
    # "                     ifelse",
    # "                 } bind dup dup",
    # "               ]",
    # "             /MatrixLMN [0.412457 0.212673 0.019334",
    # "                         0.357576 0.715152 0.119192",
    # "                         0.180437 0.072175 0.950301]",
    # "             /WhitePoint [0.9505 1.0 1.0890]",
    # "           >>",
    # "         ] setcolorspace } def",
    # "/setrgb { srgb setcolor } def",

    ### gypped code
    "/setrgb { setrgbcolor } def",
    "/rgb { setrgbcolor } def",

    "/s   { scalefont setfont } def")

    # TODO this no longer works
    # unlockBinding(".ps.prolog", asNamespace("grDevices"))
    # assignInNamespace(".ps.prolog", prolog, "grDevices")
    # lockBinding(".ps.prolog", asNamespace("grDevices"))
    cat("Opening", file, "\n")
    postscript(file, onefile, ...)
}
