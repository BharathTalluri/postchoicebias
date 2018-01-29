function MouseWait

[x, y, button] = GetMouse;
while any(button)
[x, y, button] = GetMouse;
end
while ~any(button)
[x, y, button] = GetMouse;
end