from popcorn import *
# A routine to see if two lines are touching
# [P1 P2 C1]
# y is a global

gdim = 2

PtVec = DofSpace(gdim, 0,2)
BondSca = DofSpace(1, 2,3)
GlobalVec = DofSpace(2*gdim, -1)

i_x = Input("x",PtVec)
i_y = Input("pts",GlobalVec)
i_p = Input("p",BondSca)

x1,x2 = i_x.Vertex_Split()
y1,y2 = Matrix(i_y[0:2]), Matrix(i_y[2:4])

o_test = Output("test",[BondSca],1)
s,r = symbols('s r')
itscn = solve( x1+ (x2-x1)*s - (y1+(y2-y1)*r) , s,r)
prgm = [
    IfElse( And(Ge(itscn[s],0),
                Ge(itscn[r],0),
                Le(itscn[s],1), #**2, ((x2-x1).T*(x2-x1))[0,0] ),
                Le(itscn[r],1)), #**2, ((y2-y1).T*(y2-y1))[0,0] ) ),
            Asgn(o_test, Matrix([[1]]) ),
            Asgn(o_test, Matrix([[0]]) ),
            )
    ]
Kernel("line_intersection",
       listing=prgm)

o_f = Output("R",[PtVec],1)
pforce = i_p[0] * (x2-x1)/( sqrt(((x2-x1).T*(x2-x1))[0]) )
Kernel("bond_pressure",
       listing=[
           Asgn(o_f,Matrix([pforce,-pforce]),'=')
       ])

Husk("bonds")
