extensions [nw]
turtles-own [credenceA credenceB credenceA-holding credenceB-holding preferA preferB num-success neighbor-B pulledB]
globals [num-trials groupA groupB credence-distribution subgroup-divergence subgroup-consensus-A subgroup-consensus-B overall-subgroup-divergence overall-subgroup-consensus-A overall-subgroup-consensus-B]
to setup
  clear-all
  reset-ticks
  set num-trials 1000
  create-topology
  set credence-distribution []
  ask turtles [set-credences
               set-colors]
  if count (turtles with [color = blue]) = 0
      [setup]
  ask turtles
      [ifelse count turtles with [preferA = true] >= 2
        [let stdev-A standard-deviation [credenceA] of turtles with [preferA = true]
         set overall-subgroup-consensus-A stdev-A]
        [set overall-subgroup-consensus-A 0]
       ifelse count turtles with [preferB = true] >= 2
         [let stdev-B standard-deviation [credenceB] of turtles with [preferB = true]
          set overall-subgroup-consensus-B stdev-B]
         [set overall-subgroup-consensus-B 0]
       let Aturtles [credenceA] of turtles with [preferA = true]
       let Bturtles [credenceB] of turtles with [preferB = true]
        ifelse Aturtles != [] and Bturtles != []
          [let meanA mean [credenceA] of turtles with [preferA = true]
           let meanB mean [credenceA] of turtles with [preferB = true]
           set overall-subgroup-divergence abs (meanA - meanB)]
          [if Aturtles = []
             [set overall-subgroup-divergence 0]
           if Bturtles = []
             [set overall-subgroup-divergence 0]]]
end

to create-topology
  if turtle-network = "cycle"
      [create-turtles num-turtles
       layout-circle turtles max-pxcor - 1
       ask turtles [create-links-with min-n-of 2 other turtles [distance myself]]]
  if turtle-network = "complete"
       [set-default-shape turtles "circle"
        create-turtles num-turtles
        layout-circle turtles max-pxcor - 1
        ask turtles [create-links-with other turtles]]
  if turtle-network = "preferential"
    [nw:generate-preferential-attachment turtles links num-turtles 1 [let rootturtle max-one-of turtles [count link-neighbors] layout-radial turtles links rootturtle]]
  ask turtles [set shape "dot"]
  if turtle-network = "small world"
    [nw:generate-watts-strogatz turtles links num-turtles 3 .1 repeat num-turtles [layout-tutte (turtles with [count link-neighbors > 5]) links 8]]
end

to set-credences
  set credenceA random-float 1
  if credenceA = 0 [set credenceA random-float 1]
   set credenceB (1 - credenceA)
   if credenceA > credenceB
     [set preferA true set preferB false]
   if credenceA < credenceB
     [set preferB true set preferA false]
  if credenceA = credenceB
     [set preferB false set preferA false]
end
to set-colors
  if credenceA > credenceB [set color red]
  if credenceA < credenceB [set color blue]
end

to go
  if A-better? = num-turtles [show "A is better" stop]
  if B-better? = num-turtles [show "B is better" stop]
  ask turtles [pull-the-lever-cronk]
  ask turtles [gather-others]
  ask turtles [synchronous-update]
  set groupA turtles with [preferA = true]
  set groupB turtles with [preferB = true]
  if any? turtles with [preferA = true] and any? turtles with [preferB = true] [set subgroup-consensus-A update-subgroup-consensus-A]
  if any? turtles with [preferA = true] and any? turtles with [preferB = true] [set subgroup-consensus-B update-subgroup-consensus-B]
  if any? turtles with [preferA = true] and any? turtles with [preferB = true] [set subgroup-divergence update-subgroup-divergence]
  tick
end

to-report A-better?
  report (count turtles with [credenceA > credenceB])
end

to-report B-Better?
  report count turtles with [credenceB > .99]
end


to pull-the-lever-cronk
  set pulledB false
  ifelse preferA = false and preferB = false
     [ifelse random 2 = 1
         [stop]
         [pull-B set pulledB true]]
      [ifelse preferA = true
         [stop]
         [pull-B set pulledB true]]
end

to pull-B
  set num-success 0
  repeat num-trials
    [let B random-float 1
     if B < p-success
       [set num-success 1 + num-success]]
  update-B
end

to update-B
  let ln-prior ln (credenceB / credenceA)
  let lnpEB num-success * ln (p-success) + (num-trials - num-success) * ln (1 - p-success)
  let lnpEA num-success * ln (1 - p-success) + (num-trials - num-success) * ln (p-success)
  let ln-post ln-prior + lnpEB - lnpEA
  if ln-post > 34.539 [set ln-post 34.539]
  if ln-post < -34.539 [set ln-post -34.539]
  set credenceB (1 / (1 + exp (- ln-post)))
  set CredenceA (1 / (1 + exp (ln-post)))
end

to gather-others
  if update-procedure = "Bayesian"
  [let ln-prior ln (CredenceB / CredenceA)
  let neighbor-success sum [num-success] of link-neighbors with [pulledB = true]
  let count-neighbors count link-neighbors with [pulledB = true]
  let lnpEB neighbor-success * ln (p-success) + ((num-trials * count-neighbors) - neighbor-success) * ln (1 - p-success)
  let lnpEA neighbor-success * ln (1 - p-success) + ((num-trials * count-neighbors) - neighbor-success) * ln (p-success)
  let ln-post ln-prior + lnpEB - lnpEA
  if ln-post > 34.539 [set ln-post 34.539]
  if ln-post < -34.539 [set ln-post -34.539]
  set credenceB (1 / (1 + exp (- ln-post)))
  set credenceA (1 / (1 + exp (ln-post)))]

  if update-procedure = "UPCO"
  [set neighbor-B [credenceB] of link-neighbors
  let num-neighbors length neighbor-B
  set credenceA-holding credenceA
  set credenceB-holding credenceB
  if neighbor-weights = "dependent"
  [foreach neighbor-B
     [x -> set credenceB-holding ((x ^ (1 / num-neighbors)) * credenceB-holding) / (((x ^ (1 / num-neighbors)) * credenceB-holding) + (credenceA-holding) * ((1 - x) ^ (1 / num-neighbors)))
      set credenceA-holding (1 - credenceB-holding)
      if credenceA-holding = 1 [set credenceA-holding .999999999999999 set credenceB-holding 1 - credenceA-holding]
      if credenceA-holding = 0 [set credenceA-holding .000000000000001 set credenceB-holding 1 - credenceA-holding]]]
  if neighbor-weights = "unweighted"
  [foreach neighbor-B
     [x -> set credenceB-holding (x * credenceB-holding) / ((x * credenceB-holding) + ((credenceA-holding) * (1 - x)))
      set credenceA-holding (1 - credenceB-holding)
      if credenceA-holding = 1 [set credenceA-holding .999999999999999 set credenceB-holding 1 - credenceA-holding]
      if credenceA-holding = 0 [set credenceA-holding .000000000000001 set credenceB-holding 1 - credenceA-holding]]]]
  if update-procedure = "LA"
    [set neighbor-B [credenceB] of link-neighbors
     set credenceB-holding credenceB
     set credenceB-holding ((sum neighbor-B + credenceB) / (count link-neighbors + 1))
     set credenceA-holding 1 - credenceB-holding]
end

to synchronous-update
  if update-procedure = "UPCO" or update-procedure = "LA"
      [set credenceB credenceB-holding
       set credenceA credenceA-holding]
  ifelse credenceA > credenceB
     [set color red set preferA true set preferB false]
      [ifelse credenceA < credenceB
        [set color blue set preferA false set preferB true]
        [set color green set preferA false set preferB false]]
end

to-report update-subgroup-consensus-A
  ifelse count groupA >= 2
     [let stdev-A standard-deviation [credenceA] of groupA
      let holding-subgroup-consensus-A overall-subgroup-consensus-A
      set overall-subgroup-consensus-A holding-subgroup-consensus-A + ((stdev-A - holding-subgroup-consensus-A) / (ticks + 1))
      report stdev-A]
     [report overall-subgroup-consensus-A]
end

to-report update-subgroup-consensus-B
  ifelse count groupB >= 2
     [let stdev-B standard-deviation [credenceB] of groupB
      let holding-subgroup-consensus-B overall-subgroup-consensus-B
      set overall-subgroup-consensus-B holding-subgroup-consensus-B + ((stdev-B - holding-subgroup-consensus-B) / (ticks + 1))
      report stdev-B]
     [report overall-subgroup-consensus-B]
end

to-report update-subgroup-divergence
  ifelse any? groupA and any? groupB
      [let meanA mean [credenceA] of groupA
       let meanB mean [credenceA] of groupB
       let new-subgroup-divergence abs (meanA - meanB)
       let old-subgroup-divergence overall-subgroup-divergence
       set overall-subgroup-divergence old-subgroup-divergence + ((new-subgroup-divergence - old-subgroup-divergence) / (ticks + 1))
       report abs (meanA - meanB)]
      [report overall-subgroup-divergence]
end
@#$#@#$#@
GRAPHICS-WINDOW
275
59
712
497
-1
-1
13.0
1
10
1
1
1
0
0
0
1
-16
16
-16
16
0
0
1
ticks
30.0

BUTTON
18
121
81
154
NIL
setup
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

CHOOSER
11
163
152
208
turtle-network
turtle-network
"cycle" "complete" "preferential" "small world"
2

BUTTON
89
121
152
154
NIL
go
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

PLOT
762
11
1145
186
Turtle Preferences
Ticks
Proportion of Turtles
0.0
1.0
0.0
1.0
true
true
"" ""
PENS
"Turtles that Prefer A" 1.0 0 -2674135 true "" "plot (count turtles with [preferA = true]) / num-turtles"
"Turtles that Prefer B" 1.0 0 -13345367 true "plot (count turtles with [preferB = true]) / num-turtles" "plot (count turtles with [preferB = true]) / num-turtles"

SLIDER
9
272
181
305
Num-Turtles
Num-Turtles
10
100
100.0
1
1
NIL
HORIZONTAL

CHOOSER
11
218
182
263
p-success
p-success
0.501 0.51 0.52 0.53
0

PLOT
763
201
1146
363
Subgroup Consensus
Ticks
Spread Around Group Means
0.0
10.0
0.0
1.0
true
true
"" ""
PENS
"subgroup-consensus-A" 1.0 0 -2674135 true "" "plot subgroup-consensus-A"
"subgroup-consensus-B" 1.0 0 -13345367 true "" "plot subgroup-consensus-B"

PLOT
762
378
1146
541
Subgroup Divergence
Ticks
Difference in Subgroup Means
0.0
10.0
0.0
1.0
true
true
"" ""
PENS
"subgroup-divergence" 1.0 0 -13840069 true "" "plot subgroup-divergence"

MONITOR
259
507
405
552
Divergence
overall-subgroup-divergence
17
1
11

MONITOR
417
507
565
552
Consensus A
overall-subgroup-consensus-A
17
1
11

MONITOR
575
508
725
553
Consensus B
overall-subgroup-consensus-B
17
1
11

CHOOSER
8
311
182
356
neighbor-weights
neighbor-weights
"dependent" "unweighted"
0

CHOOSER
15
69
153
114
update-procedure
update-procedure
"UPCO" "Bayesian" "LA"
0

@#$#@#$#@
## WHAT IS IT?

One procedure for pooling credences together, called UPCO, has a formal property called 'synergy.' A procedure is synergistic when, in some circumstances, the credence that results from pooling can fall outside the interval between the lowest and highest prior credence. This model investigates the conditions under which synergy, in a variety of social networks, causes polarization.

## HOW IT WORKS

Agents try to solve a bandit problem. Agents know the success rate for one arm, but do not know the success rate for the other arm. In reality, the unknown arm is better--in that pulling the arm is more likely to result in a success (1) over a failure (0)--than the known arm. Their task is to find the best arm. 

First, agents pull the arm that they're most confident is the best. If they pull the arm with the known success rate, this action is entirely uninformative, since they already know the success rate for that arm. If agents pull the arm with the unknown success rate, they do so 1000 times and privately observe the results. After collecting the results, each agent conditionalizes on their observed proportion of successes to failures drawn from the distribution for that arm. 

Next, agents pool their credences with other agents they have access to. Access is represented via undirected links between agents. Agents can pool via linear averaging or UPCO. For the sake of a comparison, agents can also perform Bayesian conditionalization on their neighbors' observations, should their neighbor pull the unknown arm. 

Agents will set their color based on which arm they think is best. If an agent is blue, they're most confident that the unknown arm is optimal. If an agent is red, they're most confident that the known arm is optimal.

This proces continues until either all agents believe the known arm is the optimal arm--in which case, no agent is ever motivated to try the unknown arm--or until all agents are extremely confident (> .99) that the unknown arm is the optimal arm. 

## HOW TO USE IT

'setup' allows you to create a network of turtles and initialize all their properties. 

'go' allows you to run the model. 

'Update-procedure' allows you to select between UPCO, linear averaging (LA), and Bayesian conditionalization for the step of the model where pooling (can) occur. 

'Turtle-network' allows you to choose between four possible social networks: cycle, complete, preferential, and small world. 

'p_success' allows you to determine the probability that pulling the unknown arm results in a success. 

'Num-turtles' determines the number of turtles in a given network, ranging from 10 to 100.

'Neighbor-weights' determines whether, in the UPCO calculation, neighbors' credences are exponentially weighted to (i) partially account for dependance, and (ii) reduce the extent to which UPCO is synergistic when pooling multiple credences. 

## WHAT THE PLOTS MEAN

'Turtle Preferences' tracks the proportion of turtles who are most confident in the known arm versus the proportion of turtles who are most confident in the unknown arm.

'Subgroup consensus' is a measure of polarization that tracks the spread of credences around a subgroup mean. Membership in a subgroup is determined by whether agents think one arm is to be preferred over another. A value closer to zero indicates very little spread around the subgroup mean, whereas a value closer to one indicates more spread around the subgroup mean. Intuitively, the less spread around subgroup means, the more polarized subgroups are. 

'Subgroup divergence' is a measure of polarization that tracks the distance between subgroup means. A value closer to zero indicates less distance between subgroup means, whereas a value closer to one indicates more distance between subgroup means. Intuitively, the more distance between subgroup means, the more polarized subgroups are. 

The monitors under the representation of the network allow you to see the values of these measures at an individual tick. 

@#$#@#$#@
default
true
0
Polygon -7500403 true true 150 5 40 250 150 205 260 250

airplane
true
0
Polygon -7500403 true true 150 0 135 15 120 60 120 105 15 165 15 195 120 180 135 240 105 270 120 285 150 270 180 285 210 270 165 240 180 180 285 195 285 165 180 105 180 60 165 15

arrow
true
0
Polygon -7500403 true true 150 0 0 150 105 150 105 293 195 293 195 150 300 150

box
false
0
Polygon -7500403 true true 150 285 285 225 285 75 150 135
Polygon -7500403 true true 150 135 15 75 150 15 285 75
Polygon -7500403 true true 15 75 15 225 150 285 150 135
Line -16777216 false 150 285 150 135
Line -16777216 false 150 135 15 75
Line -16777216 false 150 135 285 75

bug
true
0
Circle -7500403 true true 96 182 108
Circle -7500403 true true 110 127 80
Circle -7500403 true true 110 75 80
Line -7500403 true 150 100 80 30
Line -7500403 true 150 100 220 30

butterfly
true
0
Polygon -7500403 true true 150 165 209 199 225 225 225 255 195 270 165 255 150 240
Polygon -7500403 true true 150 165 89 198 75 225 75 255 105 270 135 255 150 240
Polygon -7500403 true true 139 148 100 105 55 90 25 90 10 105 10 135 25 180 40 195 85 194 139 163
Polygon -7500403 true true 162 150 200 105 245 90 275 90 290 105 290 135 275 180 260 195 215 195 162 165
Polygon -16777216 true false 150 255 135 225 120 150 135 120 150 105 165 120 180 150 165 225
Circle -16777216 true false 135 90 30
Line -16777216 false 150 105 195 60
Line -16777216 false 150 105 105 60

car
false
0
Polygon -7500403 true true 300 180 279 164 261 144 240 135 226 132 213 106 203 84 185 63 159 50 135 50 75 60 0 150 0 165 0 225 300 225 300 180
Circle -16777216 true false 180 180 90
Circle -16777216 true false 30 180 90
Polygon -16777216 true false 162 80 132 78 134 135 209 135 194 105 189 96 180 89
Circle -7500403 true true 47 195 58
Circle -7500403 true true 195 195 58

circle
false
0
Circle -7500403 true true 0 0 300

circle 2
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240

cow
false
0
Polygon -7500403 true true 200 193 197 249 179 249 177 196 166 187 140 189 93 191 78 179 72 211 49 209 48 181 37 149 25 120 25 89 45 72 103 84 179 75 198 76 252 64 272 81 293 103 285 121 255 121 242 118 224 167
Polygon -7500403 true true 73 210 86 251 62 249 48 208
Polygon -7500403 true true 25 114 16 195 9 204 23 213 25 200 39 123

cylinder
false
0
Circle -7500403 true true 0 0 300

dot
false
0
Circle -7500403 true true 90 90 120

face happy
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 255 90 239 62 213 47 191 67 179 90 203 109 218 150 225 192 218 210 203 227 181 251 194 236 217 212 240

face neutral
false
0
Circle -7500403 true true 8 7 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Rectangle -16777216 true false 60 195 240 225

face sad
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 168 90 184 62 210 47 232 67 244 90 220 109 205 150 198 192 205 210 220 227 242 251 229 236 206 212 183

fish
false
0
Polygon -1 true false 44 131 21 87 15 86 0 120 15 150 0 180 13 214 20 212 45 166
Polygon -1 true false 135 195 119 235 95 218 76 210 46 204 60 165
Polygon -1 true false 75 45 83 77 71 103 86 114 166 78 135 60
Polygon -7500403 true true 30 136 151 77 226 81 280 119 292 146 292 160 287 170 270 195 195 210 151 212 30 166
Circle -16777216 true false 215 106 30

flag
false
0
Rectangle -7500403 true true 60 15 75 300
Polygon -7500403 true true 90 150 270 90 90 30
Line -7500403 true 75 135 90 135
Line -7500403 true 75 45 90 45

flower
false
0
Polygon -10899396 true false 135 120 165 165 180 210 180 240 150 300 165 300 195 240 195 195 165 135
Circle -7500403 true true 85 132 38
Circle -7500403 true true 130 147 38
Circle -7500403 true true 192 85 38
Circle -7500403 true true 85 40 38
Circle -7500403 true true 177 40 38
Circle -7500403 true true 177 132 38
Circle -7500403 true true 70 85 38
Circle -7500403 true true 130 25 38
Circle -7500403 true true 96 51 108
Circle -16777216 true false 113 68 74
Polygon -10899396 true false 189 233 219 188 249 173 279 188 234 218
Polygon -10899396 true false 180 255 150 210 105 210 75 240 135 240

house
false
0
Rectangle -7500403 true true 45 120 255 285
Rectangle -16777216 true false 120 210 180 285
Polygon -7500403 true true 15 120 150 15 285 120
Line -16777216 false 30 120 270 120

leaf
false
0
Polygon -7500403 true true 150 210 135 195 120 210 60 210 30 195 60 180 60 165 15 135 30 120 15 105 40 104 45 90 60 90 90 105 105 120 120 120 105 60 120 60 135 30 150 15 165 30 180 60 195 60 180 120 195 120 210 105 240 90 255 90 263 104 285 105 270 120 285 135 240 165 240 180 270 195 240 210 180 210 165 195
Polygon -7500403 true true 135 195 135 240 120 255 105 255 105 285 135 285 165 240 165 195

line
true
0
Line -7500403 true 150 0 150 300

line half
true
0
Line -7500403 true 150 0 150 150

pentagon
false
0
Polygon -7500403 true true 150 15 15 120 60 285 240 285 285 120

person
false
0
Circle -7500403 true true 110 5 80
Polygon -7500403 true true 105 90 120 195 90 285 105 300 135 300 150 225 165 300 195 300 210 285 180 195 195 90
Rectangle -7500403 true true 127 79 172 94
Polygon -7500403 true true 195 90 240 150 225 180 165 105
Polygon -7500403 true true 105 90 60 150 75 180 135 105

plant
false
0
Rectangle -7500403 true true 135 90 165 300
Polygon -7500403 true true 135 255 90 210 45 195 75 255 135 285
Polygon -7500403 true true 165 255 210 210 255 195 225 255 165 285
Polygon -7500403 true true 135 180 90 135 45 120 75 180 135 210
Polygon -7500403 true true 165 180 165 210 225 180 255 120 210 135
Polygon -7500403 true true 135 105 90 60 45 45 75 105 135 135
Polygon -7500403 true true 165 105 165 135 225 105 255 45 210 60
Polygon -7500403 true true 135 90 120 45 150 15 180 45 165 90

sheep
false
15
Circle -1 true true 203 65 88
Circle -1 true true 70 65 162
Circle -1 true true 150 105 120
Polygon -7500403 true false 218 120 240 165 255 165 278 120
Circle -7500403 true false 214 72 67
Rectangle -1 true true 164 223 179 298
Polygon -1 true true 45 285 30 285 30 240 15 195 45 210
Circle -1 true true 3 83 150
Rectangle -1 true true 65 221 80 296
Polygon -1 true true 195 285 210 285 210 240 240 210 195 210
Polygon -7500403 true false 276 85 285 105 302 99 294 83
Polygon -7500403 true false 219 85 210 105 193 99 201 83

square
false
0
Rectangle -7500403 true true 30 30 270 270

square 2
false
0
Rectangle -7500403 true true 30 30 270 270
Rectangle -16777216 true false 60 60 240 240

star
false
0
Polygon -7500403 true true 151 1 185 108 298 108 207 175 242 282 151 216 59 282 94 175 3 108 116 108

target
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240
Circle -7500403 true true 60 60 180
Circle -16777216 true false 90 90 120
Circle -7500403 true true 120 120 60

tree
false
0
Circle -7500403 true true 118 3 94
Rectangle -6459832 true false 120 195 180 300
Circle -7500403 true true 65 21 108
Circle -7500403 true true 116 41 127
Circle -7500403 true true 45 90 120
Circle -7500403 true true 104 74 152

triangle
false
0
Polygon -7500403 true true 150 30 15 255 285 255

triangle 2
false
0
Polygon -7500403 true true 150 30 15 255 285 255
Polygon -16777216 true false 151 99 225 223 75 224

truck
false
0
Rectangle -7500403 true true 4 45 195 187
Polygon -7500403 true true 296 193 296 150 259 134 244 104 208 104 207 194
Rectangle -1 true false 195 60 195 105
Polygon -16777216 true false 238 112 252 141 219 141 218 112
Circle -16777216 true false 234 174 42
Rectangle -7500403 true true 181 185 214 194
Circle -16777216 true false 144 174 42
Circle -16777216 true false 24 174 42
Circle -7500403 false true 24 174 42
Circle -7500403 false true 144 174 42
Circle -7500403 false true 234 174 42

turtle
true
0
Polygon -10899396 true false 215 204 240 233 246 254 228 266 215 252 193 210
Polygon -10899396 true false 195 90 225 75 245 75 260 89 269 108 261 124 240 105 225 105 210 105
Polygon -10899396 true false 105 90 75 75 55 75 40 89 31 108 39 124 60 105 75 105 90 105
Polygon -10899396 true false 132 85 134 64 107 51 108 17 150 2 192 18 192 52 169 65 172 87
Polygon -10899396 true false 85 204 60 233 54 254 72 266 85 252 107 210
Polygon -7500403 true true 119 75 179 75 209 101 224 135 220 225 175 261 128 261 81 224 74 135 88 99

wheel
false
0
Circle -7500403 true true 3 3 294
Circle -16777216 true false 30 30 240
Line -7500403 true 150 285 150 15
Line -7500403 true 15 150 285 150
Circle -7500403 true true 120 120 60
Line -7500403 true 216 40 79 269
Line -7500403 true 40 84 269 221
Line -7500403 true 40 216 269 79
Line -7500403 true 84 40 221 269

wolf
false
0
Polygon -16777216 true false 253 133 245 131 245 133
Polygon -7500403 true true 2 194 13 197 30 191 38 193 38 205 20 226 20 257 27 265 38 266 40 260 31 253 31 230 60 206 68 198 75 209 66 228 65 243 82 261 84 268 100 267 103 261 77 239 79 231 100 207 98 196 119 201 143 202 160 195 166 210 172 213 173 238 167 251 160 248 154 265 169 264 178 247 186 240 198 260 200 271 217 271 219 262 207 258 195 230 192 198 210 184 227 164 242 144 259 145 284 151 277 141 293 140 299 134 297 127 273 119 270 105
Polygon -7500403 true true -1 195 14 180 36 166 40 153 53 140 82 131 134 133 159 126 188 115 227 108 236 102 238 98 268 86 269 92 281 87 269 103 269 113

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270
@#$#@#$#@
NetLogo 6.4.0
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
<experiments>
  <experiment name="experiment" repetitions="10000" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="30000"/>
    <metric>(count turtles with [preferB = true]) / (count turtles)</metric>
    <metric>overall-subgroup-consensus-A</metric>
    <metric>overall-subgroup-consensus-B</metric>
    <metric>overall-subgroup-divergence</metric>
    <enumeratedValueSet variable="update-procedure">
      <value value="&quot;Bayesian&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p_success">
      <value value="0.53"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-trials">
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="turtle-network">
      <value value="&quot;complete&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="neighbor-weights">
      <value value="&quot;dependent&quot;"/>
    </enumeratedValueSet>
    <steppedValueSet variable="Num-Turtles" first="10" step="10" last="50"/>
  </experiment>
  <experiment name="test" repetitions="1000" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <exitCondition>if all? turtles [credenceB &gt; .999] [show "B is better" stop]
  if all? turtles [credenceA &gt; credenceB] [show "A is better" stop]</exitCondition>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="update-procedure">
      <value value="&quot;Bayesian&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="hard-problem">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="turtle-network">
      <value value="&quot;complete&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="credence?">
      <value value="&quot;random&quot;"/>
    </enumeratedValueSet>
  </experiment>
</experiments>
@#$#@#$#@
@#$#@#$#@
default
0.0
-0.2 0 0.0 1.0
0.0 1 1.0 0.0
0.2 0 0.0 1.0
link direction
true
0
Line -7500403 true 150 150 90 180
Line -7500403 true 150 150 210 180
@#$#@#$#@
0
@#$#@#$#@
