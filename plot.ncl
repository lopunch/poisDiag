begin
    nx = fbinrecread("./p.dat", 0, 1, "integer")
    ny = fbinrecread("./p.dat", 1, 1, "integer")
    nz = fbinrecread("./p.dat", 2, 1, "integer")
;   f  = fbinrecread("./p.dat", 3, (/nz,ny,nx/), "double")
;   p  = fbinrecread("./p.dat", 4, (/nz,ny,nx/), "double")
    f  = fbinrecread("./p.dat", 3, (/nz,ny,nx/), "float")
    p  = fbinrecread("./p.dat", 4, (/nz,ny,nx/), "float")

    wks = gsn_open_wks("x11", "aaa")
    res = True
    res@gsnDraw   = False
    res@gsnFrame  = False
    res@cnFillOn  = True
    res@cnLinesOn = False
    res@cnLineLabelsOn = False
    res@cnFillMode = "CellFill"

    plot = new(2, graphic)
    resP = True
    do k = 0, nz-1, 10
        plot(0) = gsn_csm_contour(wks, f(k,:,:), res)
        plot(1) = gsn_csm_contour(wks, p(k,:,:), res)
        gsn_panel(wks, plot, (/1,2/), resP)
    end do

end
