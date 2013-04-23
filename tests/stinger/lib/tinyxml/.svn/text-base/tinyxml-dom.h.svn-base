// Copyright 2010 Christoph Gärtner
// Distributed under the Boost Software License, Version 1.0
// See <license.txt> or <http://www.boost.org/LICENSE_1_0.txt> for details

/*	an (incomplete) DOM Level 1 wrapper for tinyXML
*/

#ifndef TINYXML_DOM_H
#define TINYXML_DOM_H

typedef struct tx_node *dom_node_t;

extern size_t dom_getElementsByTagName(
	dom_node_t root, const char *name, size_t max_nodes, dom_node_t *nodes);

extern dom_node_t dom_getElementById(
	dom_node_t root, const char *id);

extern dom_node_t dom_nextSibling(
	dom_node_t node);

extern dom_node_t dom_prevSibling(
	dom_node_t node);

extern dom_node_t dom_firstChild(
	dom_node_t node);

extern _Bool dom_hasChildNodes(
	dom_node_t node);

extern size_t dom_childNodes(
	dom_node_t node, size_t max_nodes, dom_node_t *nodes);

extern const char *dom_getAttribute(
	dom_node_t node, const char *name);

#endif
