% Aufgabe 4
function x = blatt2_aufgabe4()
format long;
% Berechnung von Aufgabe 4a
a = cos(10^-8)-1

print = 'Hier beginnt Aufgabe 4b'
% Lösung 4b
% Um das Ergebnis mit einer Genauigkeit von 10^-14 berechnen zu 
% können, muss man den Cosinus mit Hilfe seiner Reihendarstellung 
% berechnen und so lange neue Glieder hinzu addieren, bis diese
% kleiner als 10^-14 sind.
x = CosinusSeries(10^-8, 14);
% ...oder z.B auf min. 50 Stellen genau
%y = CosinusSeries(10^-8, 50)
end

% Berechne den Cosinus nicht intern sondern mit seiner 
% Reihendarstellung um ein genaues Ergebnis zu erhalten
% -----------------------------------------------------
% x = Zahl deren Cosinus berechnet werden soll
% accuracy = Das Ergebnis ist auf accuracy Stellen genau
function result = CosinusSeries(x, accuracy)
i = 0;
stop = 0;
result = '';

while stop == 0
    term = (((-1)^i) * ((x^(2*i))/(factorial(2*i)))*10^accuracy);
    
    % Falls das neue Glied der Reihe zu klein, also
    % genauer ist accuracy vorgibt, kann die Rechnung
    % hier abgebrochen werden
    if (abs(term) < 1)
        stop = 1;
    end
    i = i+1;
    
    % Die einzelnen Reihenglieder werden als String
    % an 'result' angehängt. Wir rechnen die Glieder
    % nicht zusammen, weil MATLAB intern zu sehr rundet
    result = cat(2, result, ' ');
    if (term > 0)
        result = cat(2, result, '+ ');
    else 
        term = term * -1;
        result = cat(2, result, '- ');
    end
    result = cat(2, result, num2str(term*10^-accuracy));
end
end