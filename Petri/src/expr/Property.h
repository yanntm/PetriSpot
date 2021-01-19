#pragma once

#include "expr/Op.h"
#include "expr/Expression.h"


namespace petri {
namespace expr {


class Property {
	Expression * prop;
	//private PropertyType type;
	std::string name;
public :
	Property(Expression * prop, /*PropertyType type,*/ const std::string & name):
	prop(prop),name(name){
//		this.type = type;
	}
	Property() : prop(nullptr),name() {
	}
	void setName(const std::string & name) {
		this->name = name;
	}
//	public void setType(PropertyType type) {
//		this.type = type;
//	}
	void setBody(Expression * prop) {
		this->prop = prop;
	}
	const std::string & getName() {
		return name;
	}
	Expression * getBody() {
		return prop;
	}
//	public PropertyType getType() {
//		return type;
//	}
//	@Override
//	public String toString() {
//		return "Property [prop=" + prop + ", type=" + type + ", name=" + name + "]";
//	}
};



}}
