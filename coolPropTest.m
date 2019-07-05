%Cool prop test
%pyversion C:\Users\Edward\AppData\Local\Programs\Python\Python37\python.exe;
%[v,e] = pyversion; system([e,' -m pip install --user -U CoolProp']);
%Q is vapour quality
tic;
H = py.CoolProp.CoolProp.PropsSI('Hmass','P|gas',101325,'T',300,'NitrousOxide')
T = py.CoolProp.CoolProp.PropsSI('T','P|gas',101325,'Hmass',H,'NitrousOxide')
T = py.CoolProp.CoolProp.PropsSI('T','P|gas',101325,'Hmass',H,'NitrousOxide')
toc;