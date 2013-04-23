// Copyright 2010 Christoph Gärtner
// Distributed under the Boost Software License, Version 1.0
// See <license.txt> or <http://www.boost.org/LICENSE_1_0.txt> for details

/*	tinyXML is a minimal XML parser, supporting only the most basic feature set
	this excludes (among other things) document types, comments, processing
	instructions, CDATA sections and entity decoding
*/

#ifndef TINYXML_H
#define TINYXML_H

#include <stddef.h>

enum tx_types
{
	TX_EOF = 0, // only used internally
	TX_ELEMENT = 1,
	TX_ATTRIBUTE = 2,
	TX_TEXT = 4
};

struct tx_node
/*	an XML node
	element nodes have a <NULL> value
	text nodes have the name "#text"
*/
{
	enum tx_types type;
	struct tx_node *parent;
	const char *name;
	const char *value;
};

extern char *tx_parse(
	char *data, size_t max_nodes, struct tx_node *nodes);
/*	parses up to <max_nodes> of XML and returns a pointer to the unprocessed
	data; parsing is destructive, ie <data> will be modified

	returns <NULL> if <data> could be parsed completely
*/

extern struct tx_node *tx_next(
	struct tx_node *node, struct tx_node *ancestor, _Bool child,
	enum tx_types type);
/*	gets the next node in the document structure which has a given <type> and
	is a descendant of <ancestor>;  returns <NULL> if there is no such node

	if <ancestor> is <NULL>, <node->parent> will be used

	if <child> is true, only direct children of <ancestor> will be considered

	<type> may hold a combination of types by using logical or
*/


extern struct tx_node *tx_prev(
	struct tx_node *node, struct tx_node *ancestor, _Bool child,
	enum tx_types type);
/*	gets the previous node in the document structure according to the same
	semantics as <tx_next()>
*/

extern struct tx_node *tx_find(
	struct tx_node *root, struct tx_node *current, enum tx_types type,
	const char *name, const char *value, _Bool deep);
/*	finds the first child of <root> whose name or value compares equal
	to <name> or <value> respectively; returns <NULL> if no such node could
	be found

	if <current> is not <NULL>, only nodes which succeed <current> will be
	searched

	<type> may hold a combination of types by using logical or

	if either <name> or <value> is <NULL>, only the other will be considered

	if <deep> is true, all descendants of <root> (and not only direct
	children) will be considered
*/

extern struct tx_node *tx_get(
	struct tx_node *root, struct tx_node *current, size_t path_len,
	const char **path);
/*	finds a node by resolving <path>, starting at <root>; reutrns <NULL> if no
	such node could be found

	<path> must point to the first element of an array of <path_len>
	strings; each element of <path> is a node's name, where text nodes are
	represented by "#text" and attribute node names must be prefixed by "@"

	if <current> is not <NULL>, only nodes which succeed <current> will be
	considered
*/

#endif
