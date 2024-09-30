︠289327f6-17b6-44a5-a279-5573fcef62a1s︠
#The matrix: 

var('p1111','p1122','p1212','p1221','p1222','p2112','p2121','p2122','p2212','p2221','p2424','p2442','p2333','p2211','p4444','p2244','p3344','p2222','p2233','p3333','p2323','p2332')
var('q1144','q1122','q1133','q1244','q1444','q2244','q2444','q4444','q2222','q1222','q2233','q1233','q3333','q1333','q2333')

BigMatrix = matrix([[p1111      , 0          , 0          , 0          , 0          , q1144*p1122, q1122*p1122, 0          , 0          , q1133*p1122, 0          , 0          , 0          , 0          ],
                   [0          , q1144*p1212, q1144*p1221, q1244*p1222, q1244*p1222, q1444*p1222, 0          , 0          , 0          , 0          , 0          , 0          , 0          , 0          ],
                   [0          , q1144*p2112, q1144*p2121, q1244*p2122, q1244*p2122, q1444*p2122, 0          , 0          , 0          , 0          , 0          , 0          , 0          , 0          ],
                   [0          , q1244*p2212, q1244*p2221, q2244*p2424, q2244*p2442, q2444*p2333, 0          , 0          , 0          , 0          , 0          , 0          , 0          , 0          ],
                   [0          , q1244*p2212, q1244*p2221, q2244*p2442, q2244*p2424, q2444*p2333, 0          , 0          , 0          , 0          , 0          , 0          , 0          , 0          ],
                   [q1144*p2211, q1444*p2212, q1444*p2221, q2444*p2333, q2444*p2333, q4444*p4444, q2244*p2244, q1244*p2212, q1244*p2221, p3344      , 0          , 0          , 0          , 0          ],
                   [q1122*p2211, 0          , 0          , 0          , 0          , q2244*p2244, q2222*p2222, q1222*p2212, q1222*p2221, q2233*p2233, 0          , 0          , 0          , 0          ],
                   [0          , 0          , 0          , 0          , 0          , q1244*p1222, q1222*p1222, q1122*p1212, q1122*p1221, q1233*p1222, 0          , 0          , 0          , 0          ],
                   [0          , 0          , 0          , 0          , 0          , q1244*p2122, q1222*p2122, q1122*p2112, q1122*p2121, q1233*p2122, 0          , 0          , 0          , 0          ],
                   [q1133*p2211, 0          , 0          , 0          , 0          , p3344      , q2233*p2233, q1233*p2212, q1233*p2221, q3333*p3333, q1333*p2212, q1333*p2221, q2333*p2333, q2333*p2333],
                   [0          , 0          , 0          , 0          , 0          , 0          , 0          , 0          , 0          , q1333*p1222, q1133*p1212, q1133*p1221, q1233*p1222, q1233*p1222],
                   [0          , 0          , 0          , 0          , 0          , 0          , 0          , 0          , 0          , q1333*p2122, q1133*p2112, q1133*p2121, q1233*p2122, q1233*p2122],
                   [0          , 0          , 0          , 0          , 0          , 0          , 0          , 0          , 0          , q2333*p2333, q1233*p2212, q1233*p2221, q2233*p2323, q2233*p2332],
                   [0          , 0          , 0          , 0          , 0          , 0          , 0          , 0          , 0          , q2333*p2333, q1233*p2212, q1233*p2221, q2233*p2332, q2233*p2323]])
︡b8423eab-fb02-4896-9700-fd5d8328773e︡{"stdout":"(p1111, p1122, p1212, p1221, p1222, p2112, p2121, p2122, p2212, p2221, p2424, p2442, p2333, p2211, p4444, p2244, p3344, p2222, p2233, p3333, p2323, p2332)\n"}︡{"stdout":"(q1144, q1122, q1133, q1244, q1444, q2244, q2444, q4444, q2222, q1222, q2233, q1233, q3333, q1333, q2333)\n"}︡{"done":true}
︠ed12e8cd-34ed-41dc-aef1-175984057850s︠
#Find the Rank 1 Minors
Rank1MinorsMatrices = []

i=1
j=1
k=2
l=2


while i < 5:
    while k < 6:
        while j < 4: 
            while l < 5: 
                Rank1MinorsMatrices.append(BigMatrix.matrix_from_rows_and_columns([i,k], [j,l]))
                l +=1
            j +=1
            l=j+1     
        k +=1
        j=1
        l=j+1
    i += 1
    k=i+1
    j=1
    l=j+1
    
i=5
j=7
k=6
l=8


while i < 9:
    while k < 10:
        while j < 8: 
            while l < 9: 
                Rank1MinorsMatrices.append(BigMatrix.matrix_from_rows_and_columns([i,k], [j,l]))
                l +=1
            j +=1
            l=j+1     
        k +=1
        j=7
        l=j+1
    i += 1
    k=i+1
    j=7
    l=j+1

i=9
j=10
k=10
l=11


while i < 13:
    while k < 14:
        while j < 13: 
            while l < 14: 
                Rank1MinorsMatrices.append(BigMatrix.matrix_from_rows_and_columns([i,k], [j,l]))
                l +=1
            j +=1
            l=j+1     
        k +=1
        j=10
        l=j+1
    i += 1
    k=i+1
    j=10
    l=j+1

len(Rank1MinorsMatrices)
︡637bbe03-fdf0-4b2a-b72b-f1a2800a57cb︡{"stdout":"130\n"}︡{"done":true}
︠62b27322-f52e-4efe-b38e-c80dd17c9e39s︠
Rank1MinorsPolynomials = []

n=0

while n < 130:
    Rank1MinorsPolynomials.append(Rank1MinorsMatrices[n].determinant().full_simplify())
    n +=1
    
#show(Rank1MinorsPolynomials)
︡ed4a500b-d831-4fcf-a2e1-437de7830791︡{"done":true}
︠30a844d0-7249-48fc-9f75-42da50a6d08as︠
Rank2MinorsMatrices = []


j=5
k=6

while j<9:
    while k<10:
        Rank2MinorsMatrices.append(BigMatrix.matrix_from_rows_and_columns([0,j,k], [0,6,7]))
        k +=1
    j +=1
    k=j+1


i=5
j=6
k=7

while i<8:
    while j<9:
        while k<10:
            Rank2MinorsMatrices.append(BigMatrix.matrix_from_rows_and_columns([i,j,k], [0,6,7]))
            k +=1
        j +=1
        k=j+1
    i +=1
    j=i+1
    k=j+1

len(Rank2MinorsMatrices)
︡9a3a1a65-b4a2-4f71-9311-73c040bfab8b︡{"stdout":"20\n"}︡{"done":true}
︠9c83faf0-445b-4b3b-928f-75998cae21a4s︠
Rank2MinorsPolynomials = []

n=0

while n < 20:
    Rank2MinorsPolynomials.append(Rank2MinorsMatrices[n].determinant().full_simplify())
    n +=1
    
#show(Rank2MinorsPolynomials)
︡cd04105b-2c38-4c30-98dd-4faa384d8d83︡{"done":true}
︠bb28236e-d526-409d-a3ed-0b733469f9d7s︠
Rank3MinorsMatrices1 = []

j=5
k=6
l=7

while j<12:
    while k<13:
        while l<14:
            Rank3MinorsMatrices1.append(BigMatrix.matrix_from_rows_and_columns([0,j,k,l], [0,7,9,10]))
            l +=1
        k +=1
        l=k+1
    j +=1
    k=j+1
    l=k+1
    
i=5
j=6
k=7
l=8

while i<11:
    while j<12:
        while k<13:
            while l<14:
                Rank3MinorsMatrices1.append(BigMatrix.matrix_from_rows_and_columns([i,j,k,l], [0,7,9,10]))
                l +=1
            k +=1
            l=k+1
        j +=1
        k=j+1
        l=k+1
    i +=1
    j=i+1
    k=j+1
    l=k+1
    
len(Rank3MinorsMatrices1)
︡d05ae561-b6e7-49d1-8c73-b860e135154d︡{"stdout":"210\n"}︡{"done":true}
︠b2ccdb72-dbc9-4760-8685-e1b4a39951d8s︠
Rank3MinorsPolynomials1 = []

n=0

while n < 210:
    Rank3MinorsPolynomials1.append(Rank3MinorsMatrices1[n].determinant().full_simplify())
    n +=1
    
#show(Rank3MinorsPolynomials1)
︡d9ab6882-307d-4c3f-980d-eb04995ff511︡{"done":true}
︠69cfc5c1-6dcf-42a6-a8b0-2113dcfb7669s︠
Rank3MinorsMatrices2 = []

i=0
j=1
k=2
l=3

while i<7:
    while j<8:
        while k<9:
            while l<10:
                Rank3MinorsMatrices2.append(BigMatrix.matrix_from_rows_and_columns([i,j,k,l], [0,1,5,7]))
                l +=1
            k +=1
            l=k+1
        j +=1
        k=j+1
        l=k+1
    i +=1
    j=i+1
    k=j+1
    l=k+1

len(Rank3MinorsMatrices2)
︡6434eefc-aec2-404f-8440-b7b708e185b0︡{"stdout":"210\n"}︡{"done":true}
︠7d4a3548-a606-4dc5-87f1-7b6e3dad1c62s︠
Rank3MinorsPolynomials2 = []

n=0

while n < 210:
    Rank3MinorsPolynomials2.append(Rank3MinorsMatrices2[n].determinant().full_simplify())
    n +=1
    
#show(Rank3MinorsPolynomials2)
︡080e3134-f338-4fb5-bd75-bec9a9e4c642︡{"done":true}
︠d71643df-9a8f-477f-a86c-2d832b55e156s︠
var('pi1','pi2','pi3','pi4')
︡ae0c3276-4053-4489-92ef-98c1380275a3︡{"stdout":"(pi1, pi2, pi3, pi4)\n"}︡{"done":true}
︠5937a499-15f3-4f83-8c7f-afcd146342fds︠
PolysWithPis = []

m=0
while m<130:
    PolysWithPis.append(Rank1MinorsPolynomials[m].subs(q1144 = (pi1+pi2)/(pi1*pi2), q1122 = 1/((pi1+pi2)*(pi3+pi4)), q1133 = (pi3+pi4)/(pi3*pi4), q1244 = 1/(pi1*pi2), q1444 = ((pi1+pi2)*(pi2-pi1))/(pi1^2*pi2^2), q2244 = 1/(pi1*pi2*(pi3+pi4))+((pi3+pi4)-(pi1+pi2))/(pi1*pi2*(pi1+pi2)*(pi3+pi4)), q2444 = (pi2-pi1)/(pi1^2*pi2^2), q4444 = ((pi1+pi2)^2)/(pi1^2*pi2^2)+((pi1+pi2)*(pi3+pi4))/(pi1^2*pi2^2)+((pi1+pi2)*(pi1-pi2)^2)/(pi1^3*pi2^3), q2222 = 1/((pi1+pi2)^2*(pi3+pi4)^2)+(((pi1+pi2)-(pi3+pi4))^2)/((pi1+pi2)^3*(pi3+pi4)^3), q1222 = ((pi3+pi4)-(pi1+pi2))/((pi3+pi4)^2*(pi1+pi2)^2), q2233 = 1/((pi1+pi2)*pi3*pi4)+((pi1+pi2)-(pi3+pi4))/(pi3*pi4*(pi1+pi2)*(pi3+pi4)), q1233 = -1/(pi3*pi4), q3333 = ((pi3+pi4)^2)/(pi3^2*pi4^2)+((pi1+pi2)*(pi3+pi4))/(pi3^2*pi4^2)+((pi3+pi4)*(pi3-pi4)^2)/(pi3^3*pi4^3), q1333 = ((pi3+pi4)*(pi4-pi3))/(pi3^2*pi4^2), q2333 = (pi3-pi4)/(pi3^2*pi4^2)).subs(pi4=1-pi3-pi2-pi1))
    m +=1

m=0
while m<20:
    PolysWithPis.append(Rank2MinorsPolynomials[m].subs(q1144 = (pi1+pi2)/(pi1*pi2), q1122 = 1/((pi1+pi2)*(pi3+pi4)), q1133 = (pi3+pi4)/(pi3*pi4), q1244 = 1/(pi1*pi2), q1444 = ((pi1+pi2)*(pi2-pi1))/(pi1^2*pi2^2), q2244 = 1/(pi1*pi2*(pi3+pi4))+((pi3+pi4)-(pi1+pi2))/(pi1*pi2*(pi1+pi2)*(pi3+pi4)), q2444 = (pi2-pi1)/(pi1^2*pi2^2), q4444 = ((pi1+pi2)^2)/(pi1^2*pi2^2)+((pi1+pi2)*(pi3+pi4))/(pi1^2*pi2^2)+((pi1+pi2)*(pi1-pi2)^2)/(pi1^3*pi2^3), q2222 = 1/((pi1+pi2)^2*(pi3+pi4)^2)+(((pi1+pi2)-(pi3+pi4))^2)/((pi1+pi2)^3*(pi3+pi4)^3), q1222 = ((pi3+pi4)-(pi1+pi2))/((pi3+pi4)^2*(pi1+pi2)^2), q2233 = 1/((pi1+pi2)*pi3*pi4)+((pi1+pi2)-(pi3+pi4))/(pi3*pi4*(pi1+pi2)*(pi3+pi4)), q1233 = -1/(pi3*pi4), q3333 = ((pi3+pi4)^2)/(pi3^2*pi4^2)+((pi1+pi2)*(pi3+pi4))/(pi3^2*pi4^2)+((pi3+pi4)*(pi3-pi4)^2)/(pi3^3*pi4^3), q1333 = ((pi3+pi4)*(pi4-pi3))/(pi3^2*pi4^2), q2333 = (pi3-pi4)/(pi3^2*pi4^2)).subs(pi4=1-pi3-pi2-pi1))
    m +=1
    
    
m=0
while m<210:
    PolysWithPis.append(Rank3MinorsPolynomials1[m].subs(q1144 = (pi1+pi2)/(pi1*pi2), q1122 = 1/((pi1+pi2)*(pi3+pi4)), q1133 = (pi3+pi4)/(pi3*pi4), q1244 = 1/(pi1*pi2), q1444 = ((pi1+pi2)*(pi2-pi1))/(pi1^2*pi2^2), q2244 = 1/(pi1*pi2*(pi3+pi4))+((pi3+pi4)-(pi1+pi2))/(pi1*pi2*(pi1+pi2)*(pi3+pi4)), q2444 = (pi2-pi1)/(pi1^2*pi2^2), q4444 = ((pi1+pi2)^2)/(pi1^2*pi2^2)+((pi1+pi2)*(pi3+pi4))/(pi1^2*pi2^2)+((pi1+pi2)*(pi1-pi2)^2)/(pi1^3*pi2^3), q2222 = 1/((pi1+pi2)^2*(pi3+pi4)^2)+(((pi1+pi2)-(pi3+pi4))^2)/((pi1+pi2)^3*(pi3+pi4)^3), q1222 = ((pi3+pi4)-(pi1+pi2))/((pi3+pi4)^2*(pi1+pi2)^2), q2233 = 1/((pi1+pi2)*pi3*pi4)+((pi1+pi2)-(pi3+pi4))/(pi3*pi4*(pi1+pi2)*(pi3+pi4)), q1233 = -1/(pi3*pi4), q3333 = ((pi3+pi4)^2)/(pi3^2*pi4^2)+((pi1+pi2)*(pi3+pi4))/(pi3^2*pi4^2)+((pi3+pi4)*(pi3-pi4)^2)/(pi3^3*pi4^3), q1333 = ((pi3+pi4)*(pi4-pi3))/(pi3^2*pi4^2), q2333 = (pi3-pi4)/(pi3^2*pi4^2)).subs(pi4=1-pi3-pi2-pi1))
    m +=1
    
m=0
while m<210:
    PolysWithPis.append(Rank3MinorsPolynomials2[m].subs(q1144 = (pi1+pi2)/(pi1*pi2), q1122 = 1/((pi1+pi2)*(pi3+pi4)), q1133 = (pi3+pi4)/(pi3*pi4), q1244 = 1/(pi1*pi2), q1444 = ((pi1+pi2)*(pi2-pi1))/(pi1^2*pi2^2), q2244 = 1/(pi1*pi2*(pi3+pi4))+((pi3+pi4)-(pi1+pi2))/(pi1*pi2*(pi1+pi2)*(pi3+pi4)), q2444 = (pi2-pi1)/(pi1^2*pi2^2), q4444 = ((pi1+pi2)^2)/(pi1^2*pi2^2)+((pi1+pi2)*(pi3+pi4))/(pi1^2*pi2^2)+((pi1+pi2)*(pi1-pi2)^2)/(pi1^3*pi2^3), q2222 = 1/((pi1+pi2)^2*(pi3+pi4)^2)+(((pi1+pi2)-(pi3+pi4))^2)/((pi1+pi2)^3*(pi3+pi4)^3), q1222 = ((pi3+pi4)-(pi1+pi2))/((pi3+pi4)^2*(pi1+pi2)^2), q2233 = 1/((pi1+pi2)*pi3*pi4)+((pi1+pi2)-(pi3+pi4))/(pi3*pi4*(pi1+pi2)*(pi3+pi4)), q1233 = -1/(pi3*pi4), q3333 = ((pi3+pi4)^2)/(pi3^2*pi4^2)+((pi1+pi2)*(pi3+pi4))/(pi3^2*pi4^2)+((pi3+pi4)*(pi3-pi4)^2)/(pi3^3*pi4^3), q1333 = ((pi3+pi4)*(pi4-pi3))/(pi3^2*pi4^2), q2333 = (pi3-pi4)/(pi3^2*pi4^2)).subs(pi4=1-pi3-pi2-pi1))
    m +=1
︡eb23782e-6418-4371-9e07-13b8feba6290︡{"done":true}
︠a18c77d1-a8f9-44a3-8e07-b3cfd3901968s︠

len(PolysWithPis)
︡57c1c6c6-013c-4a95-b730-dd6a6d6d4eb1︡{"stdout":"570\n"}︡{"done":true}
︠c8495e3e-1a7f-4b3f-ac6a-bd467a395978s︠

S.<pi1,pi2,pi3> = PolynomialRing(QQ)
S
F = S.fraction_field()
F
R.<p1111,p1122,p1212,p1221,p1222,p2112,p2121,p2122,p2212,p2221,p2424,p2442,p2333,p2211,p4444,p2244,p3344,p2222,p2233,p3333,p2323,p2332,pi1,pi2,pi3> = PolynomialRing(F)
R
︡43088ef5-dcb0-4d4e-b095-401f09db5670︡{"stdout":"Multivariate Polynomial Ring in pi1, pi2, pi3 over Rational Field\n"}︡{"stdout":"Fraction Field of Multivariate Polynomial Ring in pi1, pi2, pi3 over Rational Field\n"}︡{"stdout":"Multivariate Polynomial Ring in p1111, p1122, p1212, p1221, p1222, p2112, p2121, p2122, p2212, p2221, p2424, p2442, p2333, p2211, p4444, p2244, p3344, p2222, p2233, p3333, p2323, p2332, pi1, pi2, pi3 over Fraction Field of Multivariate Polynomial Ring in pi1, pi2, pi3 over Rational Field\n"}︡{"done":true}
︠76a2ce64-722f-4eaf-a145-adc7e5b3e23as︠
I = R.ideal(PolysWithPis)
J = R.ideal(p1111*p1122*p1212*p1221*p1222*p2112*p2121*p2122*p2212*p2221*p2424*p2442*p2333*p2211*p4444*p2244*p3344*p2222*p2233*p3333*p2323*p2332)

︡b6a61b47-36d1-44d0-a494-4b781a301a8e︡{"done":true}
︠e73b58c7-3d33-4246-b2fb-7e1f3cd7654es︠
sat = I.saturation(J)
︡08a7d2ab-e5fd-4e90-91ba-6a01f37c191d︡{"done":true}
︠8321facc-a801-4b2e-81bb-5819e7c91e60s︠
sat[0].gens()
len(sat[0].gens())
︡f99fb054-5292-47f8-bc4d-98cee05f006a︡{"stdout":"[p2323 - p2332, (pi1^2*pi3 + 2*pi1*pi2*pi3 + pi1*pi3^2 - pi1*pi3 + pi2^2*pi3 + pi2*pi3^2 - pi2*pi3)*p2233 + (pi1^2 + 2*pi1*pi2 + 3*pi1*pi3 - 2*pi1 + pi2^2 + 3*pi2*pi3 - 2*pi2 + 3*pi3^2 - 3*pi3 + 1)*p3333 + (-pi1^2*pi3 - pi1^2 - 2*pi1*pi2*pi3 - 2*pi1*pi2 - pi1*pi3^2 - 2*pi1*pi3 + 2*pi1 - pi2^2*pi3 - pi2^2 - pi2*pi3^2 - 2*pi2*pi3 + 2*pi2 - 3*pi3^2 + 3*pi3 - 1)*p2332, (3*pi1^3*pi3 + 9*pi1^2*pi2*pi3 + 3*pi1^2*pi3^2 - 6*pi1^2*pi3 + 9*pi1*pi2^2*pi3 + 6*pi1*pi2*pi3^2 - 12*pi1*pi2*pi3 - 3*pi1*pi3^2 + 4*pi1*pi3 + 3*pi2^3*pi3 + 3*pi2^2*pi3^2 - 6*pi2^2*pi3 - 3*pi2*pi3^2 + 4*pi2*pi3 + pi3^2 - pi3)*p2222 + (pi1^3 + 3*pi1^2*pi2 + 3*pi1^2*pi3 - 2*pi1^2 + 3*pi1*pi2^2 + 6*pi1*pi2*pi3 - 4*pi1*pi2 + 3*pi1*pi3^2 - 3*pi1*pi3 + pi1 + pi2^3 + 3*pi2^2*pi3 - 2*pi2^2 + 3*pi2*pi3^2 - 3*pi2*pi3 + pi2)*p3333 + (-3*pi1^3*pi3 - pi1^3 - 9*pi1^2*pi2*pi3 - 3*pi1^2*pi2 - 3*pi1^2*pi3^2 + 3*pi1^2*pi3 + 2*pi1^2 - 9*pi1*pi2^2*pi3 - 3*pi1*pi2^2 - 6*pi1*pi2*pi3^2 + 6*pi1*pi2*pi3 + 4*pi1*pi2 - pi1*pi3 - pi1 - 3*pi2^3*pi3 - pi2^3 - 3*pi2^2*pi3^2 + 3*pi2^2*pi3 + 2*pi2^2 - pi2*pi3 - pi2 - pi3^2 + pi3)*p2332, (pi1^3*pi2*pi3^2 + 2*pi1^2*pi2^2*pi3^2 + 2*pi1^2*pi2*pi3^3 - 2*pi1^2*pi2*pi3^2 + pi1*pi2^3*pi3^2 + 2*pi1*pi2^2*pi3^3 - 2*pi1*pi2^2*pi3^2 + pi1*pi2*pi3^4 - 2*pi1*pi2*pi3^3 + pi1*pi2*pi3^2)*p3344 + (-pi1^3 - 3*pi1^2*pi2 - 3*pi1^2*pi3 + 2*pi1^2 - 3*pi1*pi2^2 - 6*pi1*pi2*pi3 + 4*pi1*pi2 - 3*pi1*pi3^2 + 3*pi1*pi3 - pi1 - pi2^3 - 3*pi2^2*pi3 + 2*pi2^2 - 3*pi2*pi3^2 + 3*pi2*pi3 - pi2)*p3333 + (pi1^3 + 3*pi1^2*pi2 + 3*pi1^2*pi3 - 2*pi1^2 + 3*pi1*pi2^2 + 6*pi1*pi2*pi3 - 4*pi1*pi2 + 3*pi1*pi3^2 - 3*pi1*pi3 + pi1 + pi2^3 + 3*pi2^2*pi3 - 2*pi2^2 + 3*pi2*pi3^2 - 3*pi2*pi3 + pi2)*p2332, (pi1^3*pi3 + 3*pi1^2*pi2*pi3 + pi1^2*pi3^2 - 3*pi1^2*pi3 + 3*pi1*pi2^2*pi3 + 2*pi1*pi2*pi3^2 - 6*pi1*pi2*pi3 - 2*pi1*pi3^2 + 3*pi1*pi3 + pi2^3*pi3 + pi2^2*pi3^2 - 3*pi2^2*pi3 - 2*pi2*pi3^2 + 3*pi2*pi3 + pi3^2 - pi3)*p2244 + (pi1^3 + 3*pi1^2*pi2 + 3*pi1^2*pi3 - 2*pi1^2 + 3*pi1*pi2^2 + 6*pi1*pi2*pi3 - 4*pi1*pi2 + 3*pi1*pi3^2 - 3*pi1*pi3 + pi1 + pi2^3 + 3*pi2^2*pi3 - 2*pi2^2 + 3*pi2*pi3^2 - 3*pi2*pi3 + pi2)*p3333 + (-pi1^3*pi3 - pi1^3 - 3*pi1^2*pi2*pi3 - 3*pi1^2*pi2 - pi1^2*pi3^2 + 2*pi1^2 - 3*pi1*pi2^2*pi3 - 3*pi1*pi2^2 - 2*pi1*pi2*pi3^2 + 4*pi1*pi2 - pi1*pi3^2 - pi1 - pi2^3*pi3 - pi2^3 - pi2^2*pi3^2 + 2*pi2^2 - pi2*pi3^2 - pi2 - pi3^2 + pi3)*p2332, (pi1^4*pi3 + pi1^3*pi2*pi3 + pi1^3*pi3^2 - 2*pi1^3*pi3 - pi1^2*pi3^2 + pi1^2*pi3 + pi1*pi2^3*pi3 + pi1*pi2*pi3^2 - pi1*pi2*pi3 + pi2^4*pi3 + pi2^3*pi3^2 - 2*pi2^3*pi3 - pi2^2*pi3^2 + pi2^2*pi3)*p4444 + (-pi1^4*pi2 - 3*pi1^3*pi2^2 - 3*pi1^3*pi2*pi3 + 2*pi1^3*pi2 - 3*pi1^2*pi2^3 - 6*pi1^2*pi2^2*pi3 + 4*pi1^2*pi2^2 - 3*pi1^2*pi2*pi3^2 + 3*pi1^2*pi2*pi3 - pi1^2*pi2 - pi1*pi2^4 - 3*pi1*pi2^3*pi3 + 2*pi1*pi2^3 - 3*pi1*pi2^2*pi3^2 + 3*pi1*pi2^2*pi3 - pi1*pi2^2)*p3333 + (pi1^4*pi2 - pi1^4*pi3 + 3*pi1^3*pi2^2 + 2*pi1^3*pi2*pi3 - 2*pi1^3*pi2 - pi1^3*pi3^2 + 2*pi1^3*pi3 + 3*pi1^2*pi2^3 + 6*pi1^2*pi2^2*pi3 - 4*pi1^2*pi2^2 + 3*pi1^2*pi2*pi3^2 - 3*pi1^2*pi2*pi3 + pi1^2*pi2 + pi1^2*pi3^2 - pi1^2*pi3 + pi1*pi2^4 + 2*pi1*pi2^3*pi3 - 2*pi1*pi2^3 + 3*pi1*pi2^2*pi3^2 - 3*pi1*pi2^2*pi3 + pi1*pi2^2 - pi1*pi2*pi3^2 + pi1*pi2*pi3 - pi2^4*pi3 - pi2^3*pi3^2 + 2*pi2^3*pi3 + pi2^2*pi3^2 - pi2^2*pi3)*p2332, p2333 - p2332, p2442 - p2332, p2424 - p2332, (pi1^2*pi3 + 2*pi1*pi2*pi3 + pi1*pi3^2 - 2*pi1*pi3 + pi2^2*pi3 + pi2*pi3^2 - 2*pi2*pi3 - pi3^2 + pi3)*p1122*p2211 + (-pi1^2 - 2*pi1*pi2 - 3*pi1*pi3 + 2*pi1 - pi2^2 - 3*pi2*pi3 + 2*pi2 - 3*pi3^2 + 3*pi3 - 1)*p1111*p3333 + (-pi1^2*pi3 + pi1^2 - 2*pi1*pi2*pi3 + 2*pi1*pi2 - pi1*pi3^2 + 5*pi1*pi3 - 2*pi1 - pi2^2*pi3 + pi2^2 - pi2*pi3^2 + 5*pi2*pi3 - 2*pi2 + 4*pi3^2 - 4*pi3 + 1)*p1111*p2332, p2122*p2221 - p2121*p2332, p1222*p2221 - p1221*p2332, p2122*p2212 - p2112*p2332, p2121*p2212 - p2112*p2221, p1222*p2212 - p1212*p2332, p1221*p2212 - p1212*p2221, p1222*p2121 - p1221*p2122, p1222*p2112 - p1212*p2122, p1221*p2112 - p1212*p2121]\n"}︡{"stdout":"19\n"}︡{"done":true}
︠b5301b35-7c86-4709-bf76-b51e47e6ad48s︠
sat[0].gen(0)
︡b2f9498f-e0af-4bf5-9a5e-46b93199805e︡{"stdout":"p2323 - p2332\n"}︡{"done":true}
︠bd1742bc-67a5-493f-9da5-2f99c70256e1s︠
sat[0].gen(1)
︡cbdc3a60-a0b7-4875-ba74-c6c3c9136aa0︡{"stdout":"(pi1^2*pi3 + 2*pi1*pi2*pi3 + pi1*pi3^2 - pi1*pi3 + pi2^2*pi3 + pi2*pi3^2 - pi2*pi3)*p2233 + (pi1^2 + 2*pi1*pi2 + 3*pi1*pi3 - 2*pi1 + pi2^2 + 3*pi2*pi3 - 2*pi2 + 3*pi3^2 - 3*pi3 + 1)*p3333 + (-pi1^2*pi3 - pi1^2 - 2*pi1*pi2*pi3 - 2*pi1*pi2 - pi1*pi3^2 - 2*pi1*pi3 + 2*pi1 - pi2^2*pi3 - pi2^2 - pi2*pi3^2 - 2*pi2*pi3 + 2*pi2 - 3*pi3^2 + 3*pi3 - 1)*p2332\n"}︡{"done":true}
︠d56841ef-d23a-4192-99ae-649c6b77a550s︠
sat[0].gen(2)
︡0438115b-dc3a-4e06-882b-d628f99d441c︡{"stdout":"(3*pi1^3*pi3 + 9*pi1^2*pi2*pi3 + 3*pi1^2*pi3^2 - 6*pi1^2*pi3 + 9*pi1*pi2^2*pi3 + 6*pi1*pi2*pi3^2 - 12*pi1*pi2*pi3 - 3*pi1*pi3^2 + 4*pi1*pi3 + 3*pi2^3*pi3 + 3*pi2^2*pi3^2 - 6*pi2^2*pi3 - 3*pi2*pi3^2 + 4*pi2*pi3 + pi3^2 - pi3)*p2222 + (pi1^3 + 3*pi1^2*pi2 + 3*pi1^2*pi3 - 2*pi1^2 + 3*pi1*pi2^2 + 6*pi1*pi2*pi3 - 4*pi1*pi2 + 3*pi1*pi3^2 - 3*pi1*pi3 + pi1 + pi2^3 + 3*pi2^2*pi3 - 2*pi2^2 + 3*pi2*pi3^2 - 3*pi2*pi3 + pi2)*p3333 + (-3*pi1^3*pi3 - pi1^3 - 9*pi1^2*pi2*pi3 - 3*pi1^2*pi2 - 3*pi1^2*pi3^2 + 3*pi1^2*pi3 + 2*pi1^2 - 9*pi1*pi2^2*pi3 - 3*pi1*pi2^2 - 6*pi1*pi2*pi3^2 + 6*pi1*pi2*pi3 + 4*pi1*pi2 - pi1*pi3 - pi1 - 3*pi2^3*pi3 - pi2^3 - 3*pi2^2*pi3^2 + 3*pi2^2*pi3 + 2*pi2^2 - pi2*pi3 - pi2 - pi3^2 + pi3)*p2332\n"}︡{"done":true}
︠45bcc08c-26f5-42e5-aeda-e65a6ad6a39fs︠
sat[0].gen(3)
︡a942b481-5920-47b9-aa38-96213ebc93c8︡{"stdout":"(pi1^3*pi2*pi3^2 + 2*pi1^2*pi2^2*pi3^2 + 2*pi1^2*pi2*pi3^3 - 2*pi1^2*pi2*pi3^2 + pi1*pi2^3*pi3^2 + 2*pi1*pi2^2*pi3^3 - 2*pi1*pi2^2*pi3^2 + pi1*pi2*pi3^4 - 2*pi1*pi2*pi3^3 + pi1*pi2*pi3^2)*p3344 + (-pi1^3 - 3*pi1^2*pi2 - 3*pi1^2*pi3 + 2*pi1^2 - 3*pi1*pi2^2 - 6*pi1*pi2*pi3 + 4*pi1*pi2 - 3*pi1*pi3^2 + 3*pi1*pi3 - pi1 - pi2^3 - 3*pi2^2*pi3 + 2*pi2^2 - 3*pi2*pi3^2 + 3*pi2*pi3 - pi2)*p3333 + (pi1^3 + 3*pi1^2*pi2 + 3*pi1^2*pi3 - 2*pi1^2 + 3*pi1*pi2^2 + 6*pi1*pi2*pi3 - 4*pi1*pi2 + 3*pi1*pi3^2 - 3*pi1*pi3 + pi1 + pi2^3 + 3*pi2^2*pi3 - 2*pi2^2 + 3*pi2*pi3^2 - 3*pi2*pi3 + pi2)*p2332\n"}︡{"done":true}
︠c9c25bd9-22d9-4d13-835d-8d4510f81d56s︠
sat[0].gen(4)
︡fe3bffec-f623-4036-b11e-fcfbb200ed0b︡{"stdout":"(pi1^3*pi3 + 3*pi1^2*pi2*pi3 + pi1^2*pi3^2 - 3*pi1^2*pi3 + 3*pi1*pi2^2*pi3 + 2*pi1*pi2*pi3^2 - 6*pi1*pi2*pi3 - 2*pi1*pi3^2 + 3*pi1*pi3 + pi2^3*pi3 + pi2^2*pi3^2 - 3*pi2^2*pi3 - 2*pi2*pi3^2 + 3*pi2*pi3 + pi3^2 - pi3)*p2244 + (pi1^3 + 3*pi1^2*pi2 + 3*pi1^2*pi3 - 2*pi1^2 + 3*pi1*pi2^2 + 6*pi1*pi2*pi3 - 4*pi1*pi2 + 3*pi1*pi3^2 - 3*pi1*pi3 + pi1 + pi2^3 + 3*pi2^2*pi3 - 2*pi2^2 + 3*pi2*pi3^2 - 3*pi2*pi3 + pi2)*p3333 + (-pi1^3*pi3 - pi1^3 - 3*pi1^2*pi2*pi3 - 3*pi1^2*pi2 - pi1^2*pi3^2 + 2*pi1^2 - 3*pi1*pi2^2*pi3 - 3*pi1*pi2^2 - 2*pi1*pi2*pi3^2 + 4*pi1*pi2 - pi1*pi3^2 - pi1 - pi2^3*pi3 - pi2^3 - pi2^2*pi3^2 + 2*pi2^2 - pi2*pi3^2 - pi2 - pi3^2 + pi3)*p2332\n"}︡{"done":true}
︠a17ddf4f-88ba-4b6b-9a3d-a2a33c7a61b5s︠
sat[0].gen(5)
︡be80e84f-e318-4cfc-88ee-a424cacddae6︡{"stdout":"(pi1^4*pi3 + pi1^3*pi2*pi3 + pi1^3*pi3^2 - 2*pi1^3*pi3 - pi1^2*pi3^2 + pi1^2*pi3 + pi1*pi2^3*pi3 + pi1*pi2*pi3^2 - pi1*pi2*pi3 + pi2^4*pi3 + pi2^3*pi3^2 - 2*pi2^3*pi3 - pi2^2*pi3^2 + pi2^2*pi3)*p4444 + (-pi1^4*pi2 - 3*pi1^3*pi2^2 - 3*pi1^3*pi2*pi3 + 2*pi1^3*pi2 - 3*pi1^2*pi2^3 - 6*pi1^2*pi2^2*pi3 + 4*pi1^2*pi2^2 - 3*pi1^2*pi2*pi3^2 + 3*pi1^2*pi2*pi3 - pi1^2*pi2 - pi1*pi2^4 - 3*pi1*pi2^3*pi3 + 2*pi1*pi2^3 - 3*pi1*pi2^2*pi3^2 + 3*pi1*pi2^2*pi3 - pi1*pi2^2)*p3333 + (pi1^4*pi2 - pi1^4*pi3 + 3*pi1^3*pi2^2 + 2*pi1^3*pi2*pi3 - 2*pi1^3*pi2 - pi1^3*pi3^2 + 2*pi1^3*pi3 + 3*pi1^2*pi2^3 + 6*pi1^2*pi2^2*pi3 - 4*pi1^2*pi2^2 + 3*pi1^2*pi2*pi3^2 - 3*pi1^2*pi2*pi3 + pi1^2*pi2 + pi1^2*pi3^2 - pi1^2*pi3 + pi1*pi2^4 + 2*pi1*pi2^3*pi3 - 2*pi1*pi2^3 + 3*pi1*pi2^2*pi3^2 - 3*pi1*pi2^2*pi3 + pi1*pi2^2 - pi1*pi2*pi3^2 + pi1*pi2*pi3 - pi2^4*pi3 - pi2^3*pi3^2 + 2*pi2^3*pi3 + pi2^2*pi3^2 - pi2^2*pi3)*p2332\n"}︡{"done":true}
︠4c442640-14ca-4cc7-8b19-0cce135f43bcs︠
sat[0].gen(6)
︡3a1c3e91-e7b4-45a0-a8a2-d4f649728ab7︡{"stdout":"p2333 - p2332\n"}︡{"done":true}
︠ef1b6cde-f536-4ef6-85bd-3f711fea1641s︠
sat[0].gen(7)
︡d656fe8e-4b4e-47b3-8c05-547c1ba0009f︡{"stdout":"p2442 - p2332\n"}︡{"done":true}
︠9ed68d1f-0db2-454c-ae5e-261684c2cc00s︠
sat[0].gen(8)
︡a820ee12-d7af-4158-9d3a-8c076ce96b00︡{"stdout":"p2424 - p2332\n"}︡{"done":true}
︠1345df85-fada-4d43-bba4-063c0d74c82es︠
sat[0].gen(9)
︡7b7b8de6-1ffc-4c66-a695-fac09e824d98︡{"stdout":"(pi1^2*pi3 + 2*pi1*pi2*pi3 + pi1*pi3^2 - 2*pi1*pi3 + pi2^2*pi3 + pi2*pi3^2 - 2*pi2*pi3 - pi3^2 + pi3)*p1122*p2211 + (-pi1^2 - 2*pi1*pi2 - 3*pi1*pi3 + 2*pi1 - pi2^2 - 3*pi2*pi3 + 2*pi2 - 3*pi3^2 + 3*pi3 - 1)*p1111*p3333 + (-pi1^2*pi3 + pi1^2 - 2*pi1*pi2*pi3 + 2*pi1*pi2 - pi1*pi3^2 + 5*pi1*pi3 - 2*pi1 - pi2^2*pi3 + pi2^2 - pi2*pi3^2 + 5*pi2*pi3 - 2*pi2 + 4*pi3^2 - 4*pi3 + 1)*p1111*p2332\n"}︡{"done":true}
︠e906795a-4267-429c-bbf9-978bfeacdacfs︠
sat[0].gen(10)
︡3a820bc1-b2bf-47f4-b1d8-9c19b89a1988︡{"stdout":"p2122*p2221 - p2121*p2332\n"}︡{"done":true}
︠71e59434-41e8-485d-925a-dd7113378a32s︠
sat[0].gen(11)
︡b7b46789-c626-4850-992b-b5aeb993fe95︡{"stdout":"p1222*p2221 - p1221*p2332\n"}︡{"done":true}
︠d6b13a52-d10a-40b5-a11d-2382c468b39fs︠
sat[0].gen(12)
︡7766079a-d3a0-4b3d-92a7-631a2f9ed7dd︡{"stdout":"p2122*p2212 - p2112*p2332\n"}︡{"done":true}
︠1b4440e0-6434-4584-8c78-7133caa860d6s︠
sat[0].gen(13)
︡944631bd-b610-42a4-b0af-4dc587912003︡{"stdout":"p2121*p2212 - p2112*p2221\n"}︡{"done":true}
︠47692bf8-45d7-47ba-9e31-1b41b08b5da2s︠
sat[0].gen(14)
︡c91c439d-7251-4969-a26c-e7f2da02af43︡{"stdout":"p1222*p2212 - p1212*p2332\n"}︡{"done":true}
︠3cdf3c5b-7cd5-4ef1-88a2-0af88fb90109s︠
sat[0].gen(15)
︡1723c8b7-b882-416a-8fbf-a197c6d51877︡{"stdout":"p1221*p2212 - p1212*p2221\n"}︡{"done":true}
︠fd0dffe7-9c99-4864-a1c7-a7519fa2adbds︠
sat[0].gen(16)
︡49bb2162-3d68-4555-ad9d-4cb94e38372c︡{"stdout":"p1222*p2121 - p1221*p2122\n"}︡{"done":true}
︠3063966d-de3f-44c9-a4df-dfc4f68e1b0ds︠
sat[0].gen(17)
︡b8ee5670-9e53-484e-a5c5-e1aa0a912de8︡{"stdout":"p1222*p2112 - p1212*p2122\n"}︡{"done":true}
︠7ec9743e-301b-4427-b983-17a5e11caddcs︠
sat[0].gen(18)
︡9c22ca84-d93b-474d-b64d-7084e1781c09︡{"stdout":"p1221*p2112 - p1212*p2121\n"}︡{"done":true}
︠b0fc6c7e-991e-41f8-b9dc-a67e6697ebf5︠









