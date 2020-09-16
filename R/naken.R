# naken.R:

# Like naken.collapse but don't collapse a vector of strings into a single string.
#
# e.g.    c("num","sqrt(num)","ord","offset(off)")
# becomes c("num","num"       "ord",        "off")

naken <- function(s)
{
    naked <- character(length(s))
    for(i in seq_along(s))
        naked[i] <- naken.collapse(s[i])
    naked
}
# Collapse s to s single string and then "naken" it
# (i.e. return only the variables in the string, separated by "+").
#
# e.g. "x1"                            becomes "x1"
#      "sqrt(x1)"                      becomes "x1"
#      "s(x1,x4,df=4)"                 becomes "x1+x4"
#      "sqrt(x1) as.numeric(x4)"       becomes "x1"
#      c("sqrt(x1)", "as.numeric(x4)") becomes "x1"
#      `x 3`                           becomes "`x 3`" (variables in backquotes unchanged)

naken.collapse <- function(s, warn.if.minus=FALSE)
{
    s <- paste.collapse(s)
    s.org <- s
    untouchable <- get.untouchable.for.naken(s)
    s <- strip.space(untouchable$s) # strip space from everything except untouchables
                                    # for "ident" gsubs below

    if(grepl("--", s, fixed=TRUE)) # '--'causes problems because '-' gets turned to '+' below
        warning0("Consecutive '-' in formula may cause problems\n         Formula:", s.org)

    # # check for "- ident" in formula (but -1 is ok)
    #
    # # commented out because this is invisible to the user, because
    # # plotmo does not plot the -ident variable
    #
    # if(warn.if.minus && grepl("\\- *[._[:alpha:]]", s)[1])
    #     warnf("plotmo will include the variable prefixed by \"-\" in the formula\n         Formula: %s", s)

    # TODO we can't ignore "-" below because of the paste0(collapse=" + ") later below
    s <- gsub("[-*/:]", "+", s)                 # replace - / * : with +

    # next two gsubs allow us to retain "x=x1" but drop "df=99" from "bs(x=x1, df=99)"

    s <- gsub("\\(._$[[:alnum:]]+=", "(", s)    # replace "(ident=" with "("
    s <- gsub("[._$[:alnum:]]+=[^,)]+", "", s)  # delete "ident=any"

    # replace ",ident" with ")+f(ident", thus "s(x0,x1)" becomes "s(x0)f(x1)"
    s <- gsub(",([._$[:alpha:]])", ")+f(\\1", s)

    regex <- "[._$[:alnum:]]*\\("
    if(grepl(regex, s)) {
        s <- gsub(regex, "", s)                 # replace ident(
        s <- gsub("[,)][^+-]*", "", s)          # remove remaining ",arg1,arg2)"
    }
    # s is now something like x1+x2, split it on "+" for further processing
    s <- strsplit(s, "+", fixed=TRUE)[[1]]

    s <- unique(s) # remove duplicates
    # remove numbers e.g. sin(x1*x2/12) is nakened to x1+x1+12, we don't want the 12
    is.num <- sapply(s, function(x) grepl("^([0-9]|\\.[0-9])", x))
    # but keep the intercept if there is one
    which1 <- which(s == "1")
    is.num[which1] <- FALSE
    s <- paste0(s[!is.num], collapse=" + ")

    replace.untouchable.for.naken(s, untouchable$replacements)
}
# In the function naken.collapse(), terms such as [string] and `string`
# must remain the same (regardless of the enclosed string).
# That is, strings in brackets or backquotes must remain untouched.
#
# This function searches for such terms, replaces them with dummies, and
# remembers where they were in  the original string (for re-replacement later).
#
# For example, if s = "x1 + x[,2] + `x 3`" we return:
#
#     out: "x1 + x!00000! + !00001!"    # note the dummies !00000! and !00001!
#
#     replacements:
#         replacement   original
#           "[00000]"     "[,2]"
#           "[00001]"    "`x 3`"

get.untouchable.for.naken <- function(s) # utility for naken
{
    # for efficiency, check for most common case (no [ or ` in s)
    if(!grepl("[\\[\`]", s)[1])
        return(list(s=s, replacements=NULL)) # no [ or ` in s

    stopifnot(length(s) == 1)

    # out and untouchables will be the returned string and table of untouchables
    # for simplicity, create untouchables as a vec and convert to a mat at the end
    out <- ""
    untouchables <- NULL

    cs <- strsplit(s, split="")[[1]] # split into individual chars for loop efficiency
    len <- length(cs)
    i <- 1
    while(i <= len) {
        c <- cs[i]
        # i==len below is for malformed strings with extra [ or ` on end
        if((c != "[" && c != "\`") || i == len) # normal character
            out <- paste0(out, c)
        else {                                  # char is [ or `, skip to matching ] or `
            istart <- i
            nestdepth <- 0
            endchar <- if(c == "[") "]" else "\`"
            for(i in (istart+1):len) {
                if(c == "[" && cs[i] == "[")
                    nestdepth <- nestdepth + 1 # nested brackets
                if(cs[i] == endchar) {
                    if(nestdepth <= 0)
                        break
                    else
                        nestdepth <- nestdepth - 1
                }
            }
            replacement <- sprint("!%05.5g!", length(untouchables) / 2)
            out <- paste0(out, replacement)
            untouchables <- c(untouchables, replacement, substr(s, istart, i))
        }
        i <- i + 1
    }
    if(length(untouchables)== 0) # malformed s="[" or s="`"
        return(list(s=s, replacements=NULL))

    replacements <- matrix(untouchables, byrow=TRUE,
                           ncol=2, nrow=length(untouchables) / 2)

    colnames(replacements) <- c("replacement", "original")

    list(s=out, replacements=replacements)
}
# undo the effect of get.untouchable.for.naken

replace.untouchable.for.naken <- function(s, replacements)
{
    for(i in seq_len(NROW(replacements)))
        s <- gsub(replacements[i, 1], replacements[i, 2], s, fixed=TRUE)
    s
}
