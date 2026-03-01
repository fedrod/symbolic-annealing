var('alpha,beta,y,z,r,s,t')
assume(alpha > 2, beta > 2)

def b(x,s=alpha,t=beta):
    c = gamma(s+t)/gamma(s)/gamma(t)
    return c*x^(s-1)*(1-x)^(t-1)
# mode
sm = solve(diff(b(x),x),x)
xmode(alpha,beta) = solve((sm[0]/x^(alpha-2)).rhs()-x,x)[0].rhs()
#show(xmode(alpha,beta))
# x0 = lhs tangent point to beta.
s0 = solve(b(x,alpha,beta)-diff(b(x,alpha,beta),x)*x,x)
#show(s0)
x0 = s0[1].rhs()
#show(x0)
# x1 = rhs tangent point to beta.
s1 = solve(b(x,alpha,beta)-diff(b(x,alpha,beta),x)*(x-1),x)
#show(s1)
x1 = Sim(s1[0].rhs())
#show(x1)
S0 = solve(x0-x,x)
S1 = solve(x1 -x,x)
#show(S0,S1)
x0 = S0[0].rhs()
x1 = S1[0].rhs()
#show(x0,x1)
X0(alpha,beta) = x0
X1(alpha,beta) = x1
alpha=20.3
beta=2.3
xm = xmode(alpha,beta)
x0 = X0(alpha,beta)
x1 = X1(alpha,beta)
r0 = b(x0,alpha,beta)/x0
r1 = b(x1,alpha,beta)/(1-x1)
# line0 meets line1 at xp
xp = solve(r0*x - r1*(1-x),x)[0].rhs().n()

p0 = plot(r0*x,x,0,xp)
p1 = plot(r1*(1-x),x,xp,1)
p2 = plot(b(x,alpha,beta),x,0,1,color="red")
#p2+p0+p1

# General Summary
var('alpha,beta')
xm = xmode(alpha,beta)
x0 = X0(alpha,beta)
x1 = X1(alpha,beta)
r0 = b(x0,alpha,beta)/x0
r1 = b(x1,alpha,beta)/(1-x1)
# line0 meets line1 at xp
xp = solve(r0*x - r1*(1-x),x)[0].rhs()
#show(x0,xp,xm,x1)
#show(xp)

####
### Example

alpha=3
beta=4.5
xm = xmode(alpha,beta)
x0 = X0(alpha,beta)
x1 = X1(alpha,beta)
r0 = b(x0,alpha,beta)/x0
r1 = b(x1,alpha,beta)/(1-x1)
# line0 meets line1 at xp
xp = solve(r0*x - r1*(1-x),x)[0].rhs()

p0 = plot(r0*x,x,0,xp)
p1 = plot(r1*(1-x),x,xp,1)
p2 = plot(b(x,alpha,beta),x,0,1,color="red")
p2+p0+p1

# General Summary
var('alpha,beta')
xm = xmode(alpha,beta)
x0 = X0(alpha,beta)
x1 = X1(alpha,beta)
r0 = b(x0,alpha,beta)/x0
r1 = b(x1,alpha,beta)/(1-x1)
# line0 meets line1 at xp
xp = solve(r0*x - r1*(1-x),x)[0].rhs()
xp
