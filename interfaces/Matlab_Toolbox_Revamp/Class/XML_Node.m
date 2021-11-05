classdef XML_Node < handle
    
    properties
        id
        time
        dt
        atol
        rtol
    end
    
    properties(Constant = true)
        lib = 'cantera_shared'
    end
    
    methods
        %% ReactorNet class constructor
        
        function x = XML_Node(name, src, wrap)
            % :param name:
            %    String name of the XML_node that should be created
            % :param src:
            %    String XML file name from which an instance of XML_Node
            %    should be created. Reads the XML tree from the input file.
            % :param wrap:
            %    Specify the ID of the XML_Node.
            % :return:
            %    Instance of class 'XML_Node'.
            
            checklib;
            
            x.id = 0;
            if nargin == 3
                x.id = wrap;
            elseif nargin == 2
                % read tree from a file
                x.id = calllib(x.lib, 'xml_get_XML_File', src, 0);
%                 if x.id < 0
%                     error(geterr);
%                 end
            elseif nargin == 1
                x.id = calllib(x.lib, 'xml_new', name);
            end
        end
        
        %% Utility methods
        
        function clear(x)
            % Clear the XML_Node object from the memory.
            checklib;
            calllib(x.lib, 'xml_del', x.id);
        end
        
        function n = addChild(root, name, val)
            % Add a child to the root.
            %
            % :param root:
            %    Instance of class 'XML_Node'.
            % :param name:
            %    String ID of the child to be added.
            % :param val:
            %    String value to be added to the child.
            % :return:
            %    Instance of class 'XML_Node'.
            
            checklib;
            calllib(root.lib, 'xml_addChild', root.id, name, val);
        end
        
        function x = build(x, file, pre)
            % Build an XML_Node in memory from an input file.
            %
            % :param x:
            %    Instance of class 'XML_Node'.
            % :param file:
            %    String input file name.
            % :param pre:
            %    Determine the method of building. If not specified or less
            %    than zero, use XML_Node::build. Otherwise, use
            %    XML_Node::get_XML_File.
            % :return:
            %    Instance of class 'XML_Node'.
            checklib;
            
            if nargin == 3 && pre > 0
                calllib(x.lib, 'xml_get_XML_file', file, 0);
            else
                calllib(x.lib, 'xml_build', x.id, file);
            end
        end
        
        function write(x, file)
            % Write XML_Node to file.
            checklib;
            calllib(x.lib, 'xml_write', x.id, file);
        end
        
        %% XML get methods
        
        function v = child(x, loc)
            % Get the child of an XML_Node.
            checklib;
            index = calllib(x.lib, 'xml_child', x.id, loc);
            v = XML_Node('', '', index);
        end
        
        function x = findbyID(root, id1)
            % Get an XML element given its ID.
            % 
            % :param root:
            %    Instance of class 'XML_Node'.
            % :param id1:
            %    String ID of the element to search for.
            % :return:
            %    Instance of class 'XML_Node'.
            checklib;
            index = calllib(root.lib, 'xml_findID', root.id, id1);
            x = XML_Node('', '', index);
        end

        function x = findbyName(root, name)
            % Get an XML element given its name.
            % 
            % :param root:
            %    Instance of class 'XML_Node'.
            % :param name:
            %    String name of the element to search for.
            % :return:
            %    Instance of class 'XML_Node'.
            checklib;
            index = calllib(root.lib, 'xml_findByName', root.id, name);
            x = XML_Node('', '', index);
        end        

        function n = nChildren(root)
            % Get the child of an XML_Node.
            checklib;
            n = calllib(root.lib, 'xml_nChildren', root.id);
        end        
        
        function a = attribute(x, key)
            % Get the XML_Node attribute with a given key. 
            %
            % :param key:
            %    String key to look up.
            % :return:
            %    Instance of class 'XML_Node'.
            
            checklib;
            a = calllib(x.lib, 'xml_attrib', x.id, key);
            % Revisit this since it returns a string! 
        end

        function v = value(x, loc)
            % Get the value at a location in an XML_Node. 
           
            checklib;
            
            if nargin == 2
                c = x.child(loc);
                index = c.id;
            else index = x.id;
            
            calllib(x.lib, 'xml_value', x.id, index);
            % Revisit this since it returns a string! 
        end        
        
    end
end

