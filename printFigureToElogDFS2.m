%%0function [] = printFigureToElogDFS2(title, image_filename, body)

    % define content
    author = 'FL2 PolariX';
    keyword = 'Experiments';
    location = 'FLASH2';
    image = base64file(image_filename);
    
    % Create elog entry as XML  
    elogXML = createElogEntryXML(author, title, 'MEASURE', keyword, body, location, image);
    
    % '-l' is needed on Mac OS to enable raw print to a generic postscript queue RK 25.11.2013
    unix(['lpr -l -Pttflog ' elogXML]);
    
%Send