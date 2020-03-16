#######################################################################################################################################################################
# Rootfinding Methods (Bisection - False Position - Fixed Point - Newton's - Aitken's Delta Squared)
# Rachel Martin
#######################################################################################################################################################################
import math
#DEFINE FUNCTION#
def tanFunc(x):                                                  #Defines a function to evaluate f(x)
    f = (math.tan((math.pi*x)) - x - 6)                          #Endpoints are (.4, .48)   --- Starting guess .48
    return f                                                     #Used in Bisection, False Position, and Aitken's

def Dfunc(x): #DERIVATIVE OF FUNCTION ABOVE#                     
    f = math.pi*((1/(math.cos(math.pi*x))**2)) -1                #Endpoints are (.4, .48)   --- Starting guess .48
    return f                                                     #Used in Newton's

def genFunc(x):                                                  #Defines a function to evaluate f(x)
    f = math.e**(-x)                                             #Starting guess recomended is 0
    return f                                                     #Used in Fixed Point

######################################################################BISECTION METHOD################################################################################
def bisection():

#GET INPUT FROM THE USER#
    valid = 1
    while (valid == 1):
        try:
            a = float(input("Enter your left endpoint: "))  #.4          #input statement for left endpoint
        except ValueError:                                               #If the input is not a number tell the user to input one
            print("Enter an integer or decimal.")                        #prints the need for an integer or float if one is not input
            continue                                                     #return to the top of the loop
        valid = 2                                                        #Ends the loop
    valid = 1                                                            #reset for new loop
    while (valid == 1):
        try:
            b = float(input("Enter your right endpoint: "))  #.48        #input statement for right endpoint
        except ValueError:                                               #If the input is not a number tell the user to input one
            print("Enter an integer or decimal.")                        #prints the need for an integer or float if one is not input
            continue
        valid = 2                                                        #Ends the loop
    valid = 1                                                            #reset for new loop
    while (valid == 1):
        try:
            Nmax = int(input("Enter max number of iterations: "))        #input statement for the max number of iterations
        except ValueError:                                               #If the input is not a number tell the user to input one
            print("Enter an integer.")                                   #prints the need for an integer if one is not input
            continue                                                     #return to the top of the loop
        valid = 2                                                        #Ends the loop


#DECLARE VARIBLES#
    E = 5*10**-5                                                         #epsilon
    z = 1                                                                #variable to copy the sign 
    sfl = math.copysign(z, tanFunc(a))                                   #records the sign of f(a)
    pa = 0.4510472588                                                    #actual root
    print(" Approximation      Error")                                     #Outputs the heading
    n = 1                                                                #Counts the number of iterations


#FOR LOOP - BISECTION#
    for i in range (Nmax):
        p = (a + b)/2                                                    #finds the midpoint/approximation
        err = abs(p - pa)                                                #calculates error
        print(n, ": %.10f    %.10f" % (p, err))                          #the approximation and the error
        if ((b-a) < 2*E):                                                #end criteria
            print("Root: %.10f" % p)                                     #prints the root
            break                                                        #breaks out of the overall loop
        sfp = math.copysign(z, tanFunc(p))                               #records the sign of f(p)
    
        if (sfl * sfp) < 0:                                              #Left side and p are opposite signs       
            b = p                                                        #sets p as right hand
        else:                                                            #p and right hand are opposite signs
            a = p                                                        #sets p as left hand
            sfl = sfp                                                    #sets sign of left to sign of p
        n += 1                                                           #adds 1 to the number of iterations
    else:
        print("Max iterations exceeded.")                                #Prints max iterations exceeded if root not found
################################################################END BISECTION METHOD##################################################################################

################################################################FALSE POSITION METHOD#################################################################################
def false():
    
#GET INPUT FROM USER#
    valid = 1
    while (valid == 1):
        try:
            a = float(input("Enter your left endpoint: "))  #.4          #input statement for left endpoint
        except ValueError:                                               #If the input is not a number tell the user to input one
            print("Enter an integer or decimal.") 
            continue
        valid = 2                                                        #Ends the loop
    valid = 1
    while (valid == 1):
        try:
            b = float(input("Enter your right endpoint: "))  #.48        #input statement for right endpoint
        except ValueError:                                               #If the input is not a number tell the user to input one
            print("Enter an integer or decimal.")
            continue
        valid = 2                                                        #Ends the loop
    valid = 1
    while (valid == 1):
        try:
            Nmax = int(input("Enter max number of iterations: "))        #input statement for the max number of iterations
        except ValueError:                                               #If the input is not a number tell the user to input one
            print("Enter an integer.")
            continue
        valid = 2
    

#DECLARE VARIABLES#
    E = 5*10**-5                                                          #convergence tolerance
    p1 = b                                                               #old p
    p2 = b                                                               #older p
    fa = tanFunc(a)                                                      #function at a
    fb = tanFunc(b)                                                      #function at b
    z = 1                                                                #variable to copy the sign
    sfa = math.copysign(z, tanFunc(a))                                   #records the sign at f(a)
    print(" Approximation    Error")                                     #Outputs the heading
    pa = 0.4510472588                                                    #actual root
    n = 1                                                                #Counts the number of iterations

####LOOP FOR FALSE POSITION####
    for i in range (Nmax):
        p = b-fb * ((b-a)/(fb-fa))                                       #approximation of the root
        err = abs(p - pa)                                                #calculates error
        print(n, ": %.10f    %.10f" % (p, err))                               #outputs the approximation and the error
        if (i > 2):                                                      #statement executes after p has been evaluated twice
            Y = (p - p1)/(p1 -p2)                                        #Lambda
            errest = abs(Y*(p-p1)/(Y-1))                                 #stopping condition
            if (errest < E):                                             #if the stopping condition is met then the root is printed
                print("Root: %.10f" % p)
                break
        fp = tanFunc(p)                                                  #evaluates the function at p
        sfp = math.copysign(z, fp)                                       #copies the sign of the function at p
        if (sfa*sfp < 0):                                                #if the left endpoint and p are opposite signs
            b = p                                                        #sets p as the right endpoint
            fb = fp                                                      #sets the right endpoint's sign as the sign of p
        else:
            a = p                                                        #sets p as the left endpoint
            fa = fp                                                      #set's the left endpoint's sign as the sign of p
        p2 = p1                                                          #saves p-1 as p-1
        p1 = p                                                           #saves p as p-1
        n += 1                                                           #adds 1 to the number of iterations
    else:
        print("Maximum number of iterations exceeded")                   #Prints max iterations exceeded if root not found
###############################################################END FALSE POSITION METHOD##############################################################################

##################################################################FIXED POINT ITERATION###############################################################################
def fixed():

##GETS INPUT FROM USER##
    valid = 1
    while (valid == 1):                                          #Checks input
        try:
            p0 = float(input("Enter initial guess: "))           #Input for initial guess     ##ALSO P-1##
        except ValueError:                                       #If the user input cannot be converted to a float then the loop will continue
            print("Enter an integer or decimal.")
            continue
        valid = 2
    valid = 1    
    while (valid == 1):
        try:                                                     #User inputs if they need linear or super linear
            Linear = float(input("If the function is super linear enter 1, otherwise enter any other integer: "))
        except ValueError:                                       #If the user input cannot be converted to a float then the loop will continue
            print("Enter 1 for super linear or any other integer for linear.")
            continue
        valid = 2
    valid = 1
    while (valid == 1):
        try:
            Nmax = int(input("Enter max number of iterations: "))        #input statement for the max number of iterations
        except ValueError:                                               #If the input is not a number tell the user to input one
            print("Enter an integer.")
            continue
        valid = 2
        
        
##DECLARE VARIBLES##
    E = 5*10**-5                                                 #stopping condition
    n = 1                                                        #counts iterations
    gp = 0                                                       #declares varible for a function of p
    pa = 0.5671432904                                            #actual root
    print(" Approximation      Error")                           #Outputs the heading
    p1 = 0.0                                                       #Initialize P-1

    for i in range (0,Nmax):
        p = genFunc(p0)                                          #evaluates the function at the initial guess
        err = abs(p - pa)
        print(n, ": %.10f    %.10f" % (p, err))                  #outputs p
        n += 1                                                   #adds to the number of iterations preformed 
        if (Linear == 1):
            if abs(p-p0) < E:                                    #SUPER LINEAR STOPPING CONDITION
                print("Root: ", round(p, 10))                    #prints the root
                break
        else:
            if (i > 2):                                          #for more than 2 iterations
                gp = ((p - p0)/(p0 - p1))
                if abs((gp/gp-1))*abs((p-p0)) < E:               #Stopping condition
                    print("Root : ", round(p, 10))               #prints root
                    break  
        p1 = p0                                                  #Pn-2
        p0 = p                                                   #Pn-1
    else:
        print("Max iterations exceeded.")                        #Prints max iterations exceeded if root not found

    
###############################################################END FIXED POINT ITERATION##############################################################################

######################################################################NEWTON'S METHOD#################################################################################
def newton():
    
####GETS INPUT FROM USER####
    valid = 1
    while (valid == 1):
        try:
            p_0 = float(input("Enter your initial guess: "))  #.48     #input statement for initial guess
        except ValueError:                                             #If the input is not a number tell the user to input one
            print("Enter an integer or decimal.")                      #Prints need for an integer or float if a letter is entered
            continue                                                   #Returns to the start of the loop
        valid = 2                                                      #Ends the loop
    valid = 1                                                          #Set to 1 for a new loop
    while (valid == 1):
        try:
            Nmax = int(input("Enter max number of iterations: "))      #input statement for the max number of iterations
        except ValueError:                                             #If the input is not a number tell the user to input one
            print("Enter an integer.")                                 #Print need for an integer if a letter or decimal is entered
            continue                                                   #Returns to start of loop
        valid = 2                                                      #Ends the loop


####DECLARE VARIABLES####
    E = 5*10**-5                                                       #epsilon
    pa = 0.4510472588                                                  #ACTUAL P
    n = 1                                                              #Counts number of iterations

    for i in range (Nmax):
        p = p_0 - ((tanFunc(p_0))/(Dfunc(p_0)))                        #evaluates function and returns a new approximation
        print(n,": %.10f" % p, end=" ")                                #outputs the approximation
        guess_actual = abs(p - pa)                                     #calculates the error
        print("    Error: %.10f" % (guess_actual))                     #outputs the approximation
        if abs(p - p_0) < E:                                           #stopping condition
            print("Root: %.10F" % p)                                   #prints the root
            break                                                      #exits loop if root is found
        p_0 = p                                                        #saves p to p-1
        n += 1                                                         #adds 1 to number of iterations
    else:
        print("Max iterations exceeded.")                              #Prints max iterations exceeded if root not found

    

####################################################################END NEWTON'S METHOD###############################################################################

#######################################################################AIKEN'S METHOD#################################################################################
def aitken():
    #GET INPUT FROM USER#
    valid = 1
    while (valid == 1):
        try:
            a = float(input("Enter your left endpoint: "))  #.4          #input statement for left endpoint
        except ValueError:                                               #If the input is not a number tell the user to input one
            print("Enter an integer or decimal.") 
            continue
        valid = 2                                                        #Ends the loop
    valid = 1
    while (valid == 1):
        try:
            b = float(input("Enter your right endpoint: "))  #.48        #input statement for right endpoint
        except ValueError:                                               #If the input is not a number tell the user to input one
            print("Enter an integer or decimal.")
            continue
        valid = 2                                                        #Ends the loop
    valid = 1
    while (valid == 1):
        try:
            Nmax = int(input("Enter max number of iterations: "))        #input statement for the max number of iterations
        except ValueError:                                               #If the input is not a number tell the user to input one
            print("Enter an integer.")
            continue
        valid = 2
    

#DECLARE VARIABLES#
    E = 5*10**-5                                                         #convergence tolerance
    p1 = 0                                                               #old p
    p2 = 0                                                               #older p
    fa = tanFunc(a)                                                      #function at a
    fb = tanFunc(b)                                                      #function at b
    z = 1                                                                #variable to copy the sign
    sfa = math.copysign(z, fa)                                           #records the sign at f(a)
    print(" Approximation    Error")                                     #Outputs the heading
    pa = 0.4510472588                                                    #actual root
    phat_1 = 0                                                           #previous value of phat
    n = 1                                                                #counts the number of iterations


###For Loop - Aitken's###
    for i in range (Nmax):
        p = b-fb * ((b-a)/(fb-fa))                                       #approximation of the root using false position
        if (i < 2):
            print(n,":")                                                 #outputs the iteration number
        if (i >= 2):                                                     #executes branch after the values needed have been calculated
            deltaPn = p -p1                                              #Delta P
            deltaSQPn = p - 2*p1 + p2                                    #Delta Squared P
            phat = p - (((deltaPn)**2)/(deltaSQPn))                      #Calculates Phat
            err = abs(phat - pa)                                         #error of the approximation and actual root
            print(n, ": %.10f    %.10f" % (phat, err))                   #outputs the approximation and the error
            if (abs(phat-phat_1) < E):                                   #Stopping condition
                print("Root: ", round(phat, 10))                         #prints the root
                break
            phat_1 = phat                                                #previous value of phat
        fp = tanFunc(p)                                                  #evaluates the function at p
        sfp = math.copysign(z, fp)                                       #copies the sign of the function at p
        if (sfa*sfp < 0):                                                #if the left endpoint and p are opposite signs
            b = p                                                        #sets p as the right endpoint
            fb = fp                                                      #sets the right endpoint's sign as the sign of p
        else:
            a = p                                                        #sets p as the left endpoint
            fa = fp                                                      #set's the left endpoint's sign as the sign of p
        p2 = p1                                                          #saves p-1 as p-2
        p1 = p                                                           #saves p as p-1
        n += 1                                                           #increases iteration number
    else:
        print("Max iterations exceeded")
####################################################################END AITKEN'S METHOD###############################################################################

###USER CHOOSES METHOD###
valid = 1
while (valid == 1):
    try:
        method = float(input("Enter 1 for Bisection Method, 2 for False Position Method, 3 for Fixed Point Iteration, 4 for Newton\'s Method, 5 for Aitken\'s Delta Squared Method: "))  #User chooses their method
        if (method > 5) or (method < 1):                                #Insures the input matches the options provided
            print("Enter a number between 1 and 4")
            continue
    except ValueError:                                                  #If the input is not a number tell the user to input one
        print("Enter an integer based on the desired method.")
        continue
    valid = 2
if (method == 1):                                                       #Calls the bisection method.
    print("You\'ve selected the Bisection Method. The recomeneded interval is (.4,.48).")
    bisection()
elif (method == 2):                                                     #Calls the false position method.
    print("You\'ve selected the Method of False Position. The recomeneded interval is (.4,.48).")
    false()
elif (method == 3):                                                     #Calls the fixed point iteration method.
    print("You\'ve selected Fixed Point Iteration. The recomeneded starting point is 0.")
    fixed()
elif (method == 4):                                                     #Calls Newton's method.
    print("You\'ve selected the Newton's Method. The recomeneded starting point .48.")
    newton()
elif (method == 5):                                                     #Calls Aitken's Delta Squared method.
    print("You\'ve selected the Aitken's Delta Squared Method. The recomeneded interval is (.4,.48).")
    aitken()