#/***********************************************************************
# * Licensed Materials - Property of IBM 
# *
# * IBM SPSS Products: Statistics Common
# *
# * (C) Copyright IBM Corp. 2014
# *
# * US Government Users Restricted Rights - Use, duplication or disclosure
# * restricted by GSA ADP Schedule Contract with IBM Corp. 
# ************************************************************************/

# author__ = "SPSS, JKP"
# version__ = "1.0.1"

# History
# 13-mar-2014 Original Version
# 11-aug-2022 check intervention variable measurement level


helptext='STATS RDD DEPENDENT=varname RUNNING=varname
    THRESHOLD=value TREATMENT=varname
    INDEP=variable list CLUSTER=varname
/OPTIONS KERNEL="triangular"* or "rectangular" or "epanechnikov"
    or "quartic" or "triweight" or "tricube" or "gaussian" or
    "cosine"
    BANDWIDTH=list of values
    SETYPE=HC1* or HC3 or const or HC or HC0 or HC1 or HC2 or HC4
    or HC4m or HC5
    MCCRARY=YES* or NO
    PLOT = YES* or NO

DEPENDENT, RUNNING, and THRESHOLD are required. Values marked with *
are defaults.

DEPENDENT names the dependent or outcome variable.
RUNNING names the running or assignment variable.
THRESHOLD specifies the cutoff above which the treatment is given.

If TREATMENT is not specified, the design is sharp; if a variable
is specified, the design is fuzzy.

INDEP is a list of covariates.

CLUSTER specifies a variable whose values are cluster numbers within
which the errors are assumed correlated.  If specified, robust
standard errors are calculated regardless of the SETYPE setting.
The R package documentation suggests that if the running variable
is discrete that the data be clustered by its values.

KERNEL specifies the kernel to be used.

BANDWIDTH specifies a list of bandwidths.  If not specified,
the bandwidth is calculated by the Imbens-Kalyanaraman method, and
that is used along with half and double the value.

SETYPE specifies the robust standard error method to be used if
not overridden by the use of CLUSTER.  CONST assumes constant
variances while the others assume heteroscedasticity.  HC0
is White\'s estimator.  The others are variations of it.

MCCRARY=YES produces the McCrary sorting test

PLOT=YES produces a plot of the outcome vs the running variable and,
for fuzzy designs, a plot of the treatment vs the running variable.
The plot setting also controls whether the McCrary test produces
a plot.

STATS RDD /HELP.  prints this information and does nothing else.
'

### MAIN ROUTINE ###
dordd = function(dep, running, threshold, indep=NULL, treatment=NULL,
    cluster=NULL, kernel="triangular", bandwidth=NULL, setype="hc1", mccrary=TRUE, doplot=TRUE
    ) {
    #estimate regression discontinuity model
    
    setuplocalization("STATS_RDD")
    
    # A warnings proc name is associated with the regular output
    # (and the same omsid), because warnings/errors may appear in
    # a separate procedure block following the regular output
    procname=gtxt("Regression Discontinuity Analysis")
    warningsprocname = gtxt("Regression Discontinuity: Warnings")
    omsid="STATSRDD"
    warns = Warn(procname=warningsprocname,omsid=omsid)
    if (setype != "const") {
        setype = toupper(setype)
    }
    if (setype == "HC4M") {
        setype = "HC4m"
    }
    allargs = as.list(environment())

    tryCatch(library(rdd), error=function(e){
        warns$warn(gtxtf("The R %s package is required but could not be loaded.", "rdd"),dostop=TRUE)
        }
    )
    
    frml = buildfrml(dep, running, indep, treatment, warns)
    # need to expose these variables to plot for later :-( 
    assign("frml",frml, envir=.GlobalEnv)
    assign("threshold",threshold, envir=.GlobalEnv)
    assign("cluster",cluster, envir=.GlobalEnv)
    assign("kernel",kernel, envir=.GlobalEnv)
    assign("bandwidth", bandwidth, envir=.GlobalEnv)
    assign("setype",setype, envir=.GlobalEnv)
    
    alldata = c(dep, running, indep, treatment, cluster)
    dta = spssdata.GetDataFromSPSS(alldata, missingValueToNA=TRUE, factorMode="levels")
    if (any(sapply(dta[c(2:(2+length(indep)))], is.factor))) {
      print('oops')
      warns$warn("The intervention variable and covariates must have a scale measurement level", dostop=TRUE)

    }

    if (is.null(cluster)) {
        clustersetting = NULL
    } else {
        clustersetting = dta[[cluster]]
    }
    assign("clustersetting", clustersetting, envir=.GlobalEnv)
    res = tryCatch(RDestimate(formula=frml, data=dta, cutpoint=threshold, kernel=kernel, 
            bw=bandwidth, se.type=setype, cluster=clustersetting),
        error = function(e) {
            warns$warn(e$message, dostop=TRUE)
            return(NULL)
        }
    )

    # RDsummary uses cat.  Need to suppress that output
    f = tempfile()
    sink(f)
    ressum = tryCatch(summary(res), error=function(e) {return(NULL)})
    sink()
    unlink(f)
    
    displayresults(allargs, dta, res, ressum, warns)
    
    # wrap up, flushing warnings and clearing workspace
    warns$warn()  # ensure warnings are displayed
    rm(list=ls(envir=.GlobalEnv))
    rm(list=ls())
}

buildfrml = function(dep, running, indep, treatment, warns) {
    # Return formula expression as formula object

    # dep is the name of dependent variable
    # running is the name of the running (assignment) variable
    # indep is the list of names of the independent variables
    # treatment is the name of the treatment variable
    # warns is the error message object
    
   if (!is.null(indep)) {
       cov = paste(indep, collapse="+")
       cov = paste("|", cov)
   } else {
       cov = ""
   }
   if (is.null(treatment)) {
       frml = paste(dep, "~", running, cov, collapse=" ")
   } else {
       frml = paste(dep, "~", paste(running, treatment, sep ="+"), cov, collapse=" ")
   }
   return(as.formula(frml))
}


displayresults = function(allargs, dta, res, ressum, warns) {
    # Produce pivot tables and charts
    
    StartProcedure(allargs[["procname"]], allargs[["omsid"]])
    summarylabels=list(
        gtxt("Dependent Variable"),
        gtxt("Running Variable"),
        gtxt("Threshold"),
        gtxt("Model Type"),
        gtxt("Endogenous Treatment Variable"),
        gtxt("Kernel"),
        gtxt("Std. Error Type"),
        gtxt("Cluster Variable")
    )
    
    # NOTE: some terms below are from syntax and would need to mapped for
    # translation purposes, e.g., the censortype and ties parameters
    summaryvalues = list(
        allargs[["dep"]],
        allargs[["running"]],
        allargs[["threshold"]],
        res$type,
        ifelse(is.null(allargs[["treatment"]]), gtxt("--None--"), allargs[["treatment"]]),
        allargs[["kernel"]],
        ifelse(is.null(allargs[["cluster"]]), allargs[["setype"]], gtxt("Robust")),
        ifelse(is.null(allargs[["cluster"]]), gtxt("--None--"), allargs[["cluster"]])
    )

    names(summaryvalues) = summarylabels
    summarydf = data.frame(cbind(summaryvalues))
    colnames(summarydf) = gtxt("Values")

    spsspivottable.Display(summarydf, title=gtxt("Regression Discontinuity Summary"), 
                           templateName="STATSRDDSUMMARY",
                           caption=gtxt("Results computed by R rdd package"),
                           isSplit=FALSE,
                           format=formatSpec.Count
    )

    if (!"coefficients" %in% names(ressum)) {
        warns$warn(gtxt("Model could not be estimated"))
    } else {
        # coefficients
        coefs = data.frame(ressum$coefficients)
        names(coefs) = c(gtxt("Bandwidth"), 
            gtxt("Number of Cases"), 
            ifelse(res$type == "sharp", gtxt("Average Treatment\nEffect Estimate"), gtxt("Wald Estimate")),
            gtxt("Std. Error"), gtxt("Z Value"), gtxt("Sig."))
        spsspivottable.Display(coefs, title=gtxt("Coefficients"),
            template="STATSRDDCOEF",
            isSplit=FALSE
        )

        # F statistics
        fstat = data.frame(ressum$fstat)
        names(fstat) = c(gtxt("F"), gtxt("Numerator D.F."), gtxt("Denom. D.F."), gtxt("Sig."))
        spsspivottable.Display(fstat, title=gtxt("F Statistics"), 
            template="STATSRDDF",
            isSplit=FALSE
        )
        
        if (allargs[["doplot"]]) {
            plot(res, which=1)
            title(main=gtxt("Outcome Variable Vs. Running Variable"),
                  xlab=allargs[["running"]], ylab=allargs[["dep"]])
            if (res$type == "fuzzy") {
                plot(res, which=2)
                title(main=gtxt("Treatment Variable Vs. Running Variable"),
                      xlab=allargs[["running"]], ylab=allargs[["treatment"]])
            }
        }
        
        if (allargs[["mccrary"]]) {
            mcc = DCdensity(runvar=dta[[allargs[["running"]]]], 
                cutpoint=allargs[["threshold"]], ext.out=TRUE, plot=allargs[["doplot"]])
            if (allargs[["doplot"]]) {
                title(main=gtxt("McCrary Sorting Test"), 
                    sub=sprintf(gtxt("Running Variable: %s"), allargs[["running"]]))
            }
            lbls = list(
                sprintf(gtxt("Estimated Log Difference at Threshold (%s)"), round(allargs[["threshold"]], 3)),
                gtxt("Std. Error"),
                gtxt("Z Statistic"),
                gtxt("Sig."),
                gtxt("Bin Size"),
                gtxt("Bandwidth")
            )
    
            values = list(
                round(mcc$theta, 5),
                round(mcc$se, 5),
                round(mcc$z, 5),
                round(mcc$p, 5),
                round(mcc$binsize, 5),
                round(mcc$bw, 5)
            )
            names(values) = lbls
            mccdf = data.frame(cbind(values))
            colnames(mccdf) = gtxt("Statistics")
            spsspivottable.Display(mccdf, title=gtxt("McCrary Sorting Test"),
                template="STATSRDDMCCRARY",
                caption = sprintf(gtxt("Running Variable: %s"), allargs[["running"]]),
                isSplit=FALSE)
        }
    }
    spsspkg.EndProcedure()
}



Warn = function(procname, omsid) {
    # constructor (sort of) for message management
    lcl = list(
        procname=procname,
        omsid=omsid,
        msglist = list(),  # accumulate messages
        msgnum = 0
    )
    # This line is the key to this approach
    lcl = list2env(lcl) # makes this list into an environment

    lcl$warn = function(msg=NULL, dostop=FALSE, inproc=FALSE) {
        # Accumulate messages and, if dostop or no message, display all
        # messages and end procedure state
        # If dostop, issue a stop.

        if (!is.null(msg)) { # accumulate message
            assign("msgnum", lcl$msgnum + 1, envir=lcl)
            # There seems to be no way to update an object, only replace it
            m = lcl$msglist
            m[[lcl$msgnum]] = msg
            assign("msglist", m, envir=lcl)
        } 

        if (is.null(msg) || dostop) {
            lcl$display(inproc)  # display messages and end procedure state
            if (dostop) {
                stop(gtxt("End of procedure"), call.=FALSE)  # may result in dangling error text
            }
        }
    }
    
    lcl$display = function(inproc=FALSE) {
        # display any accumulated messages as a warnings table or as prints
        # and end procedure state, if any

        if (lcl$msgnum == 0) {   # nothing to display
            if (inproc) {
                spss.EndProcedure()
            }
        } else {
            if (!inproc) {
                procok =tryCatch({
                    StartProcedure(lcl$procname, lcl$omsid)
                    TRUE
                    },
                    error = function(e) {
                        FALSE
                    }
                )
            }
            if (procok) {  # build and display a Warnings table if we can
                table = spss.BasePivotTable("Warnings ","Warnings") # do not translate this
                rowdim = BasePivotTable.Append(table,Dimension.Place.row, 
                    gtxt("Message Number"), hideName = FALSE,hideLabels = FALSE)

                for (i in 1:lcl$msgnum) {
                    rowcategory = spss.CellText.String(as.character(i))
                    BasePivotTable.SetCategories(table,rowdim,rowcategory)
                    BasePivotTable.SetCellValue(table,rowcategory, 
                        spss.CellText.String(lcl$msglist[[i]]))
                }
                spsspkg.EndProcedure()   # implies display
            } else { # can't produce a table
                for (i in 1:lcl$msgnum) {
                    print(lcl$msglist[[i]])
                }
            }
        }
    }
    return(lcl)
}

# localization initialization
setuplocalization = function(domain) {
    # find and bind translation file names
    # domain is the root name of the extension command .R file, e.g., "SPSSINC_BREUSCH_PAGAN"
    # This would be bound to root location/SPSSINC_BREUSCH_PAGAN/lang

    fpath = Find(file.exists, file.path(.libPaths(), paste(domain, ".R", sep="")))
    bindtextdomain(domain, file.path(dirname(fpath), domain, "lang"))
} 
# override for api to account for extra parameter in V19 and beyond
StartProcedure <- function(procname, omsid) {
    if (substr(spsspkg.GetSPSSVersion(),1, 2) >= 19) {
        spsspkg.StartProcedure(procname, omsid)
    }
    else {
        spsspkg.StartProcedure(omsid)
    }
}

gtxt <- function(...) {
    return(gettext(...,domain="STATS_RDD"))
}

gtxtf <- function(...) {
    return(gettextf(...,domain="STATS_RDD"))
}


Run = function(args) {
    #Execute the STATS COXREGR command

    cmdname = args[[1]]
    args = args[[2]]
    oobj = spsspkg.Syntax(list(
        spsspkg.Template("DEPENDENT", subc="",  ktype="existingvarlist", var="dep",
            islist=FALSE),
        spsspkg.Template("RUNNING", subc="", ktype="existingvarlist", var="running",
            islist=FALSE),
        spsspkg.Template("THRESHOLD", subc="",  ktype="float", var="threshold"),
        spsspkg.Template("INDEP", subc="", ktype="existingvarlist", var="indep", islist=TRUE),
        spsspkg.Template("TREATMENT", subc="", ktype="existingvarlist", var="treatment"),
        spsspkg.Template("CLUSTER", subc="", ktype="existingvarlist", var="cluster", islist=FALSE),
        
        spsspkg.Template("KERNEL", subc="OPTIONS", ktype="str", var="kernel",
            vallist=list("triangular", "rectangular", "epanechnikov", "quartic", "triweight",
            "tricube", "gaussian", "cosine")),
        spsspkg.Template("BANDWIDTH", subc="OPTIONS", ktype="float", var="bandwidth", islist=TRUE,
            vallist=list(0)),
        spsspkg.Template("SETYPE", subc="OPTIONS", ktype="str", var="setype",
            vallist=list("hc1", "const", "hc3", "hc0", "hc", "hc2", "hc4", "hc4m", "hc5")),
        spsspkg.Template("MCCRARY", subc="OPTIONS", ktype="bool", var="mccrary"),
        spsspkg.Template("PLOT", subc="OPTIONS", ktype="bool", var="doplot")
    ))

    # A HELP subcommand overrides all else
    if ("HELP" %in% attr(args,"names")) {
        #writeLines(helptext)
        helper(cmdname)
    }
    else {
        res <- spsspkg.processcmd(oobj, args, "dordd")
    }
}

helper = function(cmdname) {
    # find the html help file and display in the default browser
    # cmdname may have blanks that need to be converted to _ to match the file
    
    fn = gsub(" ", "_", cmdname, fixed=TRUE)
    thefile = Find(file.exists, file.path(.libPaths(), fn, "markdown.html"))
    if (is.null(thefile)) {
        print("Help file not found")
    } else {
        browseURL(paste("file://", thefile, sep=""))
    }
}
if (exists("spsspkg.helper")) {
assign("helper", spsspkg.helper)
}
