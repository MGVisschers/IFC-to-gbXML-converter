# Overview
This IFC to gbXML converter is a result of research performed during my graduation project at the Eindhoven University of Technology. The converter aims to create an improved interoperability between BIM and whole-building energy analysis by translating the open exchange format IFC to a validated gbXML file format. Relationships between gbXML elements and corresponding IFC entities are created and linkages are established by iterating over multiple IFC entities.
#
Author: Maarten Visschers
#
# Using the converter
The current version of the converter is Python based without a GUI, with the script you are able to run the IFC-to-gbXML conversion. Supported are:  
* The IFC2X3_TC1 schema (http://www.buildingsmart-tech.org/ifc/IFC2x3/TC1/html/index.htm);  
* The gbXML 6.01 schema (http://www.gbxml.org/Schema_Current_GreenBuildingXML_gbXML). 

#
**What to do:**   
For the operation of the converter it is needed to include 2nd level space boundaries (*IfcRelSpaceBoundary*) in the IFC file. This entity is used to create explicit geometry for the gbXML schema. More details can be found by using the provided test cases and studying the thesis.  
#
**Capabilities:** 
* Be compatible with IFC2x3_TC1 files;  
* Process basic IFC building elements (e.g. walls and floors);  
* Convert implicit to explicit geometric information (to comply with gbXML);  
* Process IFC window and opening elements;  
* Write valid information according to the official gbXML schema;  
* Create gbXML files which are suitable for whole-building energy simulation in DesignBuilder;  
* Include gbXML construction types with each their corresponding layer and materials;  
* Include thermal properties (e.g. IFC analytical property sets).  

# More information
For more information, please contact maartenvisschers@hotmail.com, or file an issue in the repository.
