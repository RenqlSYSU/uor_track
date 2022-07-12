begin
    a="abc"
    print(a)
    print([/"b","c"/])
    if(a.eq."a") then
        print(a)
    end if

    if(a.eq.[/"bcef","cdfrg"/]) then
        print(a+"in")
    end if
end
