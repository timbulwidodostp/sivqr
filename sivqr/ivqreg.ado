*! ivqreg version DowonKwak 1.0.0 28July2010
program define ivqreg, eclass byable(onecall)
        version 10.1

        if _by() {
                local BY `"by `_byvars'`_byrc0':"'
        }
        if replay() {
                if ("`e(cmd)'"!="myprog") error 301

                Replay `0'
                exit
        }
        `BY' Estimate `0'
end

program Estimate, eclass byable(recall) sort

local cmdline `"ivqreg `0'"'

* step 1: parsing
gettoken lhs 0 : 0

gettoken p 0 : 0, parse(" (") quotes
while `"`p'"' != "(" {
	local inexog `"`inexog' `p'"'
	gettoken p 0 : 0, parse(" (") quotes
}

gettoken q 0 : 0, parse(" =") quotes
while `"`q'"' != "=" {
	local endo `"`endo' `q'"'
	gettoken q 0 : 0, parse(" =") quotes
}

gettoken r 0 : 0, parse(" )") quotes
while `"`r'"' != ")" {
	local instr `"`instr' `r'"'
	gettoken r 0 : 0, parse(" )") quotes
}

local allinstr `"`inexog' `instr'"'
local rhs `"`endo' `inexog'"'

* step 2: Syntax

local options "Level(cilevel)"
        syntax [if] [in] [,                ///
        `options'                          ///
         noConstant                        /// 
         Quantil(real 0.5)                 ///
         Robust                            /// 
         * ]
 		if "`constan'"!="" {
                        di in red "nocons invalid"
                        exit 198
                }
                if `quantil' >= 1 {
                        local quant = `quantil'/100
                }
                else    local quant "`quantil'"
                if `quant' <= 0 | `quant' >= 1 {
                        di in red "quantiles(`quantil') out of range"
                        exit 198
                }
	
	/* count what is left over after marking out sample             */
        cap count if `touse'
        local N = r(N)
        if (`N' == 0) error 2000

* Step 3: Initial Two stage quantile regression for construnction of grid
marksample touse
markout `touse' `lhs' `endo' `instr' `inexog' `allinstr' `rhs' 

qui regress `endo' `allinstr' if `touse'
predict dhat, xb
local dhat dhat
keep if `touse'

display _newline "Initial Estimation: `quantil'th Two Stage Quantile Regression" _col(60) "Number of obs = " e(N)
qui qreg `lhs' `dhat' `inexog' if `touse', q(`quantil')
predict iner if `touse', resid
matrix btwo=e(b)
scalar bdht=btwo[1,1]
matrix vtwo=e(V)
local vnames `endo' `inexog' _cons
matrix rownames vtwo = `vnames'
matrix colnames vtwo = `vnames'
matrix colnames btwo = `vnames'
ereturn post btwo vtwo, depname(`lhs') obs(`e(N)') esample(`touse')
ereturn display

display _newline in yellow "Grid search is in progress (200)"

marksample touse
markout `touse'
tempname intv nn nz dh_se

mata: ivqr_initial("`lhs'", "`endo'", "`inexog'", "`instr'", "`dhat'", "iner", "`touse'")
drop iner
local nz nz 
local nn nn

*set of grid
scal dh_se=sqrt((`nn'/100)*intv[1,1])
scal gll=bdht-2*dh_se
scal size=dh_se/50
matrix gr_alpha=J(200,2,.)
	
forvalues j=1/200 {
matrix gr_alpha[`j',1]=gll+`j'*size
}  

* step 3: Grid search		
forvalues j=1/50 {
	gen yo`j'k=`lhs'-`endo'*gr_alpha[`j',1]
	qui qreg yo`j'k `instr' `inexog' if `touse', q(`quantil')
	matrix b_al`j'l=e(b)
	predict er`j', resid
	display _col(`j') _continue .
}
display _col(52) 50

forvalues j=51/100 {
	gen yo`j'k=`lhs'-`endo'*gr_alpha[`j',1]
	qui qreg yo`j'k `instr' `inexog' if `touse', q(`quantil')
	matrix b_al`j'l=e(b)
	predict er`j', resid
	local k= `j'-50
	display _col(`k') _continue .
}

display _col(52) 100


forvalues j=101/150 {
	gen yo`j'k=`lhs'-`endo'*gr_alpha[`j',1]
	qui qreg yo`j'k `instr' `inexog' if `touse', q(`quantil')
	matrix b_al`j'l=e(b)
	predict er`j', resid
	local k= `j'-100
	display _col(`k') _continue .
}

display _col(52) 150

forvalues j=151/200 {
	gen yo`j'k=`lhs'-`endo'*gr_alpha[`j',1]
	qui qreg yo`j'k `instr' `inexog' if `touse', q(`quantil')
	matrix b_al`j'l=e(b)
	predict er`j', resid
	local k= `j'-150
	display _col(`k') _continue .
}

display _col(52) 200


* step 4: IV quantile regression 		
forvalues j=1/200 {
marksample touse
markout `touse'
tempname vq
matrix b_gam=b_al`j'l[1,1..`nz']
mata: ivqr_grid("yo`j'k", "`endo'", "`inexog'", "`instr'", "er`j'", "`touse'")
matrix gr_alpha[`j',2] = vq
drop yo`j'k er`j'
		  }

	svmat gr_alpha, name(a_grid)	
	gsort a_grid2
 	scal alpha_opt=a_grid1	  
	drop a_grid1 a_grid2

* step 5: Standard errors for just-identified case
if `nz'==1 {      
   gen double yoptm = `lhs'-`endo'*alpha_opt
   qui qreg  yoptm  `instr' `inexog' if `touse', q(`quantil')
   matrix xzhat=e(b)
   matrix xhat=xzhat[.,`nz'+1...]
   qui predict ero, resid
   qui iqreg ero
   matrix iqrr=e(b)
   qui sum ero
   scal eeh = r(sd)		
   scal iqq = iqrr[1,1]/1.349
   scal iqr = min(iqq,eeh)
   matrix bhat=(alpha_opt, xhat)                                                        
   *This is the final estimate for theta
   matrix colnames bhat = `endo' `inexog' _cons

marksample touse
markout `touse'
tempname Vrb vn
mata: ivqr_rbstse("yoptm", "`endo'", "`inexog'", "`instr'", "ero", "`touse'")
drop ero dhat yoptm	
}	

* step 5-2: Standard errors for over-identified case
if `nz'>1 {

   gen double yoptm = `lhs'-`endo'*alpha_opt
   qui qreg  yoptm  `instr' `inexog' if `touse', q(`quantil')
   matrix xzhat=e(b)
   matrix xhat=xzhat[.,`nz'+1...]
   qui predict ero, resid
   qui iqreg ero
   matrix iqrr=e(b)
   qui sum ero
   scal eeh = r(sd)		
   scal iqq = iqrr[1,1]/1.349
   scal iqr = min(iqq,eeh)			
   matrix bhat=(alpha_opt, xhat)                                                        
   *This is the final estimate for theta
   matrix colnames bhat = `endo' `inexog' _cons


marksample touse
markout `touse' 
tempname Vrb vn
mata: ivqrovi_rbstse("yoptm", "`endo'", "`inexog'", "`instr'", "ero", "`touse'")

drop ero dhat yoptm
}	

* step 6: return saved results
*// row and column names
local vnames `endo' `inexog' _cons
matrix rownames vn = `vnames'
matrix colnames vn = `vnames'
matrix rownames Vrb = `vnames'
matrix colnames Vrb = `vnames'
matrix colnames bhat = `vnames'

*// robust standard error option
local N = r(N)
if "`robust'" !="" {
	ereturn post bhat Vrb, depname(`lhs') obs(`N') esample(`touse')
	
		   }
	else    ereturn post bhat vn, depname(`lhs') obs(`N') esample(`touse')
		 
		   
*// Return as e() macro
ereturn local depvar ="`lhs'"
ereturn scalar N = r(N)
ereturn local cmdline `cmdline'
ereturn local title "Instrumental Variable `quantil'th Quantile Regression"
ereturn local cmd "ivqreg"

display _newline "`quantil'th Instrumental Variable Quantile Regression" _col(60) "Number of obs = " e(N)
ereturn display

end
