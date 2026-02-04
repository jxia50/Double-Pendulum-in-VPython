from vpython import *

## intialize graphs
energyGraph = graph(title = "Energy vs Time", xtitle = "Time", ytitle = "Energy", width = 600, height = 300, align = 'right')
peCurve = gcurve(graph = energyGraph, color = color.red, label = "PE", fast = True)
keCurve = gcurve(graph = energyGraph, color = color.blue, label = "KE", fast = True)
teCurve = gcurve(graph = energyGraph, color = color.green, label = "Total Energy", fast = True)

## initialize parameters
m1 = 1.0 ## mass of joint 1
m2 = 1.0 ## mass of joint 2
l1 = 1.0 ## length of rod 1
l2 = 1.0 ## length of rod 2
g = 9.81 ## graivty

theta1 = pi ## angle of rod 1 from vertical
theta2 = pi/6 ## angle of rod 2 from vertical
omega1 = 0.0 ## angular velocity of joint 1
omega2 = 0.0 ## angular velocity of joint 2

dt = 0.0001 ## time step
t = 0.0 ## current time
scene.width = 800
scene.height = 600
scene.background = color.white
scene.title = "Double Pendulum Simulation using RK4 solver"
scene.align = 'left'
origin = vector(0, 0, 0)
frameRate = 3000
trails = True
scene.userzoom = False

## initialize objects
pivot = sphere(pos = origin, radius = 0.05, color = color.red)
rod1 = cylinder(pos = origin, radius = 0.02, color = color.blue)
joint1 = sphere(radius = 0.05 * sqrt(m1), color = color.red, make_trail = True)
rod2 = cylinder(radius = 0.02, color = color.blue)
joint2 = sphere(radius = 0.05 * sqrt(m2), color = color.red, make_trail = True)

def getPositionsFromAngles(theta1, theta2): ## get the position of joints based on angles
    x1 = l1 * sin(theta1)
    x2 = x1 + l2 * sin(theta2)
    y1 = -l1 * cos(theta1)
    y2 = y1 - l2 * cos(theta2)
    return vector(x1, y1, 0), vector(x2, y2, 0)

def updatePositions(): ## update the position of rods and joints based off current angles
    firstJointPos = getPositionsFromAngles(theta1, theta2)[0]
    secondJointPos = getPositionsFromAngles(theta1, theta2)[1]

    x1 = firstJointPos.x 
    x2 = secondJointPos.x
    y1 = firstJointPos.y
    y2 = secondJointPos.y ## where 1 is joint 1 and 2 is joint 2

    rod1.axis = vector(x1, y1, 0)
    joint1.pos = rod1.pos + rod1.axis
    rod2.pos = joint1.pos
    rod2.axis = vector(x2 - x1, y2 - y1, 0)
    joint2.pos = rod2.pos + rod2.axis

def enableTrail(button):
    global trails
    trails = not trails
    joint1.make_trail = trails
    joint2.make_trail = trails
    joint1.clear_trail()
    joint2.clear_trail()
    button.text = "Disable Trails" if trails else "Enable Trails"
    button.background = vector(1, 0.734 ,0.740) if trails else vector(0.699, 1.000,0.740)

def sliderUpdates(event): ## update parameters based on slider input
    global m1, m2, l1, l2, theta1, theta2, omega1, omega2, frameRate
    if event.id == 'm1':
        m1 = event.value
        joint1.radius = 0.05 * sqrt(m1)
    elif event.id == 'm2':
        m2 = event.value
        joint2.radius = 0.05 * sqrt(m2)
    elif event.id == 'l1':
        l1 = event.value
    elif event.id == 'l2':
        l2 = event.value
    elif event.id == 'rod1':
        theta1 = event.value
    elif event.id == 'rod2':
        theta2 = event.value
    elif event.id == 'speed':
        frameRate = event.value
        joint1.clear_trail()
        joint2.clear_trail()
        return
    
    ## reset angular velocities and trails to clean up the visuals
    joint1.clear_trail()
    joint2.clear_trail()
    omega1 = 0.0
    omega2 = 0.0

    ## reset energy graphs
    peCurve.delete()
    keCurve.delete()
    teCurve.delete()

def PE(): ## get the potential energy of the system
    y1 = -l1 * cos(theta1)
    y2 = y1 - l2 * cos(theta2)
    U = m1 * g * y1 + m2 * g * y2
    return U

def KE(): ## get the kinetic energy of the system
    v1_sq = (l1 * omega1)**2
    v2_sq = (l1 * omega1)**2 + (l2 * omega2)**2 + 2 * l1 * l2 * omega1 * omega2 * cos(theta1 - theta2)
    K = 0.5 * m1 * v1_sq + 0.5 * m2 * v2_sq 
    return K

def updateGraph(): ## update the graphs of various energy vs time (pe, ke, total vs. time)
    global t, dt
    pe = PE()
    ke = KE()

    if(int(t * (1 / dt)) % 250 == 0): ## only update graph every 250 time steps
        peCurve.plot(pos=(t, pe))
        keCurve.plot(pos=(t, ke))
        teCurve.plot(pos=(t, pe + ke))

def lagrangeRhs(theta1, theta2, omega1, omega2): ## get the right hand side of the lagrange equations using formulas for angular acceleration
    delta = theta2 - theta1 ## difference in theta angles
    sin_d = sin(delta) ## precompute sin and cos of delta
    cos_d = cos(delta)

    denom = m1 + m2*sin_d*sin_d ## get the denominator term used for formulas below

    dth1 = omega1 ## derivative of theta1 is omega1
    dth2 = omega2 ## derivative of theta2 is omega2

    ## NOTE: the follow formulas are based off a report from University of Maryland which will be listed in the report
    num1 = (m2 * l1 * omega1 * omega1 * sin_d * cos_d ## numerator for domega1
            + m2 * g * sin(theta2) * cos_d
            + m2 * l2 * omega2 * omega2 * sin_d
            - (m1 + m2) * g * sin(theta1))
    domega1 = num1 / (l1 * denom)

    num2 = (-m2 * l2 * omega2 * omega2 * sin_d * cos_d ## numerator for domega2
            + (m1 + m2) * (g*sin(theta1) * cos_d - l1 * omega1 * omega1 * sin_d - g * sin(theta2)))
    domega2 = num2 / (l2 * denom)

    return [dth1, dth2, domega1, domega2] ## return derivatives of theta1, theta2, omega1, omega2

## OLD IMPLEMENTATION BELOW (for reference); this implementation makes assumptions about equal lengths and masses

    # alpha = m2 / m1
    # beta = l2 / l1
    # gamma = g / l1

    # a1 = (beta) * (m2 / (m1 + m2)) * cos(theta1 - theta2)
    # a2 = (1 / beta) * cos(theta1 - theta2)

    # f1 = -(beta) * (m2 / (m1 + m2)) * omega2**2 * sin(theta1 - theta2) - gamma * sin(theta1)
    # f2 = (1 / beta) * omega1**2 * sin(theta1 - theta2) - gamma * sin(theta2)

    # g1 = (f1 - a1 * f2) / (1 - a1 * a2)
    # g2 = (f2 - a2 * f1) / (1 - a1 * a2)

    # return [omega1, omega2, g1, g2]

def iterateTime(): ## RK4 integrator to update angles and angular velocities
    global theta1, theta2, omega1, omega2, dt
    originalValues = [theta1, theta2, omega1, omega2]
    nextValuesOfInterest = originalValues.copy()

    for i in range(4): ## loop for k1, k2, k3, k4 (intermediate slope estimates used in RK4--fourth order Runge-Kutta method)
        currentK = lagrangeRhs(nextValuesOfInterest[0], nextValuesOfInterest[1], nextValuesOfInterest[2], nextValuesOfInterest[3]) ## calculate the next k slope based on the prior k slope
        ## loop for each variable (theta1, theta2, omega1, omega2) to create 
        # the next set of values to be used to determine the next k slope, and to build R
        for j in range(4): 
            if i == 0:
                if j == 0: ## Only need to initialize R once (for i == 0 and j == 0)
                    ## R is used to solve Ordinary Differential Equations (ODEs) 
                    # by taking weighted averages of slopes at different 
                    # points (start, middle, end) within a time step
                    R = [0, 0, 0, 0]
                R[j] = 1 / 6 * dt * currentK[j] ## 1/6 * dt of k1
            if i < 2:
                if i == 1:
                    R[j] += 1 / 3 * dt * currentK[j] ## 1/3 * dt of k2 (and k3)
                nextValuesOfInterest[j] = originalValues[j] + dt * currentK[j] / 2
            elif i == 2:
                R[j] += 1 / 3 * dt * currentK[j] ## 1/3 * dt of k3 (and k2)
                nextValuesOfInterest[j] = originalValues[j] + dt * currentK[j] 
            else:
                R[j] += 1 / 6 * dt * currentK[j] ## 1/6 * dt of k4
    
    ## update theta1, theta2, omega1, omega2 using R
    theta1 += R[0]
    theta2 += R[1]
    omega1 += R[2]
    omega2 += R[3]

def main():
    updatePositions() ## initial update to positions

    ## create the sliders for user input
    ## change rate of program
    wtext(text = "Adjust program rate:\n\n")
    wtext(text = 'Program rate\n')
    speedSlider = slider(bind = sliderUpdates, min = 3000, max = 30000, id = 'speed', align = 'left', value = frameRate)
    t0 = wtext(text = 'Frame Rate = '+'{:1.2f}x'.format(speedSlider.value / 3000))
    wtext(text = '\n')
    button(text="Disable Trails", pos=scene.caption_anchor, bind=enableTrail, background = vector(1, 0.734 ,0.740))
    ## first block of sliders for starting angles
    wtext(text = '\n\n\n')
    wtext(text = "Adjust starting angle (CCW with respect to the vertical):\n\n")
    wtext(text = 'Starting angle for rod 1\n')
    angleSlider1 = slider(bind = sliderUpdates, min = 0, max = 2 * pi, id = 'rod1', align = 'left', value = theta1)
    t1 = wtext(text = 'θ1 = '+'{:1.2f}'.format(angleSlider1.value))
    wtext(text = '\nStarting angle for rod 2\n')
    angleSlider2 = slider(bind = sliderUpdates, min = 0, max = 2 * pi, id = 'rod2', align = 'left', value = theta2)
    t2 = wtext(text = 'θ2 = '+'{:1.2f}'.format(angleSlider2.value))
    wtext(text = '\n\n\n\n')

    ## second block of sliders for mass and length
    wtext(text = "Adjust parameters:\n\n")
    wtext(text = 'Mass of joint 1\n')
    mass1Slider = slider(bind = sliderUpdates, min = 0.25, max = 10, id = 'm1', align = 'left', value = m1)
    t3 = wtext(text = 'm1 = '+'{:1.2f}'.format(mass1Slider.value) + 'kg')
    wtext(text = '\nMass of joint 2\n')
    mass2Slider = slider(bind = sliderUpdates, min = 0.25, max = 10, id = 'm2', align = 'left', value = m2)
    t4 = wtext(text = 'm2 = '+'{:1.2f}'.format(mass2Slider.value) + 'kg')
    wtext(text = '\nLength of Rod 1\n')
    length1Slider = slider(bind = sliderUpdates, min = 0.25, max = 3, id = 'l1', align = 'left', value = l1)
    t5 = wtext(text = 'l1 = '+'{:1.2f}'.format(length1Slider.value) + 'm')
    wtext(text = '\nLength of Rod 2\n')    
    length2Slider = slider(bind = sliderUpdates, min = 0.25, max = 3, id = 'l2', align = 'left', value = l2)
    t6 = wtext(text = 'l2 = '+'{:1.2f}'.format(length2Slider.value) + 'm')
    wtext(text = '\n\n')


    ## main loop
    while True:
        global t
        rate(frameRate) ## set a rate for the simulation
        
        ## update time, positions, and graphs
        iterateTime()
        updatePositions()
        updateGraph()

        ## update text displays for current slider values
        t0.text = 'Frame Rate = '+'{:1.2f}x'.format(speedSlider.value / 3000)
        t1.text = 'θ1 = '+'{:1.2f}'.format(angleSlider1.value)
        t2.text = 'θ2 = '+'{:1.2f}'.format(angleSlider2.value)
        t3.text = 'm1 = '+'{:1.2f}'.format(mass1Slider.value) + 'kg'
        t4.text = 'm2 = '+'{:1.2f}'.format(mass2Slider.value) + 'kg'
        t5.text = 'l1 = '+'{:1.2f}'.format(length1Slider.value) + 'm'
        t6.text = 'l2 = '+'{:1.2f}'.format(length2Slider.value) + 'm'
    
        t = t + dt
        pass

main()
