def strChange(string,n = 6):
    tol = 1e-12
    string = string.replace(" ","")
    function_list = string.split(",")
    func = function_list[len(function_list)-1]
    funclist = []
    for i in range(len(function_list)-1):
        valueArgument = function_list[i].split("=")
        try:
            if abs(float(valueArgument[1])) < tol:
                valueArgument[1] = "0"
        except:
            pass
        for j in range(len(funclist)):
            valueArgument2 = funclist[j]
            valueArgument[1] = valueArgument[1].replace(valueArgument2[0],valueArgument2[1])
        funclist.append(valueArgument)
    for i in range(len(funclist)):
        valueArgument = funclist[i]
        func = func.replace(valueArgument[0],valueArgument[1])
    for i in range(n):
        SIN = "sin(q{})".format(i+1)
        COS = "cos(q{})".format(i+1)
        sn = "s{}".format(i+1)
        cn = "c{}".format(i+1)
        func = func.replace(SIN,sn)
        func = func.replace(COS,cn)
    return func        
